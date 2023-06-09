# get_scaling_factors.jl
#=
This code will read in XY data generated by the compare_xy_sac.jl code
as text files and calculate scaling factors in small windows between:
1) TC+Q330 vs TC+Drum (to assess recorder response)
2) TC+Q330 vs PE+Q330 (to assess PE response)
3) TC+Q330 vs PE+Drum (to assess both PE and recorder response)
test) TC+Q330 vs TC+Q330 (should be 1)
This data will then be plotted against amplitude to get an idea of the
amplitude response / instrument response of the recorder or PE.
This is all based on an assumption that TC+Q330 is "correct"

Written By: Thomas Lee
Created On: Jan. 24, 2023

Last Modified: Feb. 9, 2023
    - Added sgram plot for scl values

=#

## PACKAGES
push!(LOAD_PATH, "/Users/thomaslee/Research/Julia/MyModules")
import LeeFunctions
const lf = LeeFunctions
import Dates
using Plots
using Seis
using Interpolations
using StatsBase
using ProgressBars
using DSP
using Measures

## SETUP
# point to timed xy data
TC_DRUM_in = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220713/DigitSeis_20220713_TCZ_ASL_timed.txt"
PE_DRUM_in = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_PEZ_20220713/DigitSeis_20220713_PEZ_ASL_timed.txt"
# TC_DRUM_in = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220714/DigitSeis_20220714_TCZ_ASL_timed.txt"
# PE_DRUM_in = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_PEZ_20220714/DigitSeis_20220714_PEZ_ASL_timed.txt"
# point to sac data
TC_Q330_in = "/Users/thomaslee/Research/Microform/ASL_NewRecords/ASL_SAC/00_TCZ/"
PE_Q330_in = "/Users/thomaslee/Research/Microform/ASL_NewRecords/ASL_SAC/10_PEZ/"
# first-stage filtering options
TC_Q330_Filt = true # filter the sac to emulate an analog filter
PE_Q330_Filt = false
b1 = 0.01 # lowband
b2 = 1.34 # highband
designmethod = DSP.Butterworth(2)
# bandpass filtering options
Band_Filt_All = true
global b1_all = 0.05 # lowband 20s
global b2_all = 0.1 # highband 10s
global designmethod_all = DSP.Butterworth(2)
# windowing options
sml_wndw_size = Dates.Second(180) # must be even
sml_wndw_step = Dates.Second(30)
# culling option
lower_percentile = 0.01 # amplitude culling 
upper_percentile = 0.97 # 1 = 100%
# small plot options
plot_small_segs = true
small_seg_wndw = Dates.Minute(60)
small_seg_step = Dates.Minute(15)
# seismogram plot options
plot_sgrams = true # make a plot of seismograms with timing annotated
# output options
data_out = "/Users/thomaslee/Desktop/RESULTS_scaling_factors/"
suffix = "13"

## SETUP OUTPUT
if !isdir(data_out)
    mkdir(data_out)
end

## READ DATA
# read sac
print(string("Reading TC+Q330 sac from ",TC_Q330_in,"\n"))
# check if is dir or sac
Stmp = lf.combsac(TC_Q330_in)
TC_Q330_name = string(Stmp.sta.sta,".",Stmp.sta.cha)
TC_Q330_dt = round(Stmp.delta,digits=2)
TC_Q330_T = lf.gettime(Stmp)
TC_Q330_D = trace(Stmp)
TC_Q330_D = TC_Q330_D .- mean(filter(!isnan,TC_Q330_D)) # demean simple
Stmp = [] # clear Stmp
print(string("  Read in ",length(TC_Q330_D)," samples\n"))
print("\n")
print(string("Reading PE+Q330 sac from ",PE_Q330_in,"\n"))
# check if is dir or sac
Stmp = lf.combsac(PE_Q330_in)
PE_Q330_name = string(Stmp.sta.sta,".",Stmp.sta.cha)
PE_Q330_dt = round(Stmp.delta,digits=2)
PE_Q330_T = lf.gettime(Stmp)
PE_Q330_D = trace(Stmp)
PE_Q330_D = PE_Q330_D .- mean(filter(!isnan,PE_Q330_D)) # demean simple
Stmp = [] # clear Stmp
print(string("  Read in ",length(PE_Q330_D)," samples\n"))
print("\n")
#filter if needed
if TC_Q330_Filt
    print(string("Filtering the TC+Q330 data between ",b1," and ",b2,"Hz\n\n"))
    # get good indices
    local gidx = findall(.!isnan.(TC_Q330_D))
    # trim ends if they are NaN
    TC_Q330_D = TC_Q330_D[gidx[1]:gidx[end]]
    TC_Q330_T = TC_Q330_T[gidx[1]:gidx[end]]
    # get the gidx again after trimming
    gidx = findall(.!isnan.(TC_Q330_D))
    # get NaN
    nan_idx = findall(isnan.(TC_Q330_D))
    if !isempty(nan_idx)
    # interpolate
        no_nan_itp = LinearInterpolation(
            Dates.value.(TC_Q330_T[gidx].-TC_Q330_T[1]),
            TC_Q330_D[gidx]
        )
        TC_Q330_D = no_nan_itp(Dates.value.(TC_Q330_T.-TC_Q330_T[1]))
    end
    # filter
    tmpfilt = DSP.digitalfilter(
        DSP.Bandpass(b1, b2; fs = 1/TC_Q330_dt),
        designmethod,
    )
    TC_Q330_D = DSP.filt(
        tmpfilt, 
        TC_Q330_D,
    )
    # reintroduce NaN
    if !isempty(nan_idx)
        TC_Q330_D[nan_idx] .= NaN
    end
    # plot filter response
    w = range(0,(0.5/TC_Q330_dt),length=250) # get w vector to Nyquist since response will be normalized
    fresp = DSP.freqz(tmpfilt) # frequency response
    presp = DSP.phasez(tmpfilt) # phase response
    hpresp = plot(
        w,
        abs.(fresp),
        title = "Filter Response",
        label = "Freq",
        ylabel = "Freq Response",
        lc = :blue,
    )
    plot!(
        twinx(hpresp),
        w,
        presp,
        label = "Phase",
        ylabel ="Phase Response",
        lc = :red,
        legend = :bottomright,
    )
    plot!(
        hpresp,
        right_margin = 10mm,
    )
    savefig(hpresp,string(data_out,TC_Q330_name,"_response_",suffix,".pdf"))
end
if PE_Q330_Filt
    print(string("Filtering the PE+Q330 data between ",b1," and ",b2,"Hz\n\n"))
    # get good indices
    local gidx = findall(.!isnan.(PE_Q330_D))
    # trim ends if they are NaN
    PE_Q330_D = PE_Q330_D[gidx[1]:gidx[end]]
    PE_Q330_T = PE_Q330_T[gidx[1]:gidx[end]]
    # get tthe gidx again after trimming
    gidx = findall(.!isnan.(PE_Q330_D))
    # get NaN
    nan_idx = findall(isnan.(PE_Q330_D))
    if !isempty(nan_idx)
    # interpolate
        no_nan_itp = LinearInterpolation(
            Dates.value.(PE_Q330_T[gidx].-PE_Q330_T[1]),
            PE_Q330_D[gidx]
        )
        PE_Q330_D = no_nan_itp(Dates.value.(PE_Q330_T.-PE_Q330_T[1]))
    end
    # filter
    tmpfilt = DSP.digitalfilter(
        DSP.Bandpass(b1, b2; fs = 1/PE_Q330_dt),
        designmethod,
    )
    PE_Q330_D = DSP.filt(
        tmpfilt, 
        PE_Q330_D,
    )
    # reintroduce NaN
    if !isempty(nan_idx)
        PE_Q330_D[nan_idx] .= NaN
    end
    # plot filter response
    w = range(0,(0.5/PE_Q330_dt),length=250) # get w vector to Nyquist since response will be normalized
    fresp = DSP.freqz(tmpfilt) # frequency response
    presp = DSP.phasez(tmpfilt) # phase response
    hpresp = plot(
        w,
        abs.(fresp),
        title = "Filter Response",
        label = "Freq",
        ylabel = "Freq Response",
        lc = :blue,
    )
    plot!(
        twinx(hpresp),
        w,
        presp,
        label = "Phase",
        ylabel ="Phase Response",
        lc = :red,
        legend = :bottomright,
    )
    plot!(
        hpresp,
        right_margin = 10mm,
    )
    savefig(hpresp,string(data_out,PE_Q330_name,"_response_",suffix,".pdf"))
end
# read xy for PE_DRUM and TC_DRUM
print(string("Reading PE+DRUM xy from ",PE_DRUM_in,"\n"))
f = open(PE_DRUM_in)
# read file
ln = readlines(f)
# close file
close(f)
# setup variables
PE_DRUM_D = Vector{Float64}(undef,length(ln)-4)
PE_DRUM_T = Vector{Dates.Millisecond}(undef,length(ln)-4)
PE_DRUM_X = Vector{Float64}(undef,length(ln)-4)
PE_DRUM_Y = Vector{Float64}(undef,length(ln)-4)
PE_DRUM_Z = Vector{Float64}(undef,length(ln)-4)
PE_DRUM_stime = Dates.DateTime(ln[1][9:end])
# skip 5 headers and loop
for i = 5:lastindex(ln)
    # get commas
    commas = findall(",",ln[i])
    # set values
    PE_DRUM_D[i-4] = parse(Float64,ln[i][commas[4][end]+1:commas[5][1]-1])
    PE_DRUM_T[i-4] = Dates.Millisecond(ln[i][commas[1][end]+1:commas[2][1]-1])
    PE_DRUM_X[i-4] = parse(Float64,ln[i][commas[5][end]+1:commas[6][1]-1])
    PE_DRUM_Y[i-4] = parse(Float64,ln[i][commas[6][end]+1:commas[7][1]-1])
    PE_DRUM_Z[i-4] = parse(Float64,ln[i][commas[7][end]+1:commas[8][1]-1])
end
PE_DRUM_T = PE_DRUM_T .+ PE_DRUM_stime
print(string("  Read in ",length(PE_DRUM_D)," samples\n"))
# same for TC
print(string("Reading TC+DRUM xy from ",PE_DRUM_in,"\n"))
f = open(TC_DRUM_in)
# read file
ln = readlines(f)
# close file
close(f)
# setup variables
TC_DRUM_D = Vector{Float64}(undef,length(ln)-4)
TC_DRUM_T = Vector{Dates.Millisecond}(undef,length(ln)-4)
TC_DRUM_X = Vector{Float64}(undef,length(ln)-4)
TC_DRUM_Y = Vector{Float64}(undef,length(ln)-4)
TC_DRUM_Z = Vector{Float64}(undef,length(ln)-4)
TC_DRUM_stime = Dates.DateTime(ln[1][9:end])
# skip 5 headers and loop
for i = 5:lastindex(ln)
    # get commas
    commas = findall(",",ln[i])
    # set values
    TC_DRUM_D[i-4] = parse(Float64,ln[i][commas[4][end]+1:commas[5][1]-1])
    TC_DRUM_T[i-4] = Dates.Millisecond(ln[i][commas[1][end]+1:commas[2][1]-1])
    TC_DRUM_X[i-4] = parse(Float64,ln[i][commas[5][end]+1:commas[6][1]-1])
    TC_DRUM_Y[i-4] = parse(Float64,ln[i][commas[6][end]+1:commas[7][1]-1])
    TC_DRUM_Z[i-4] = parse(Float64,ln[i][commas[7][end]+1:commas[8][1]-1])
end
TC_DRUM_T = TC_DRUM_T .+ TC_DRUM_stime
print(string("  Read in ",length(TC_DRUM_D)," samples\n"))

## MAKE SEISMOGRAM PLOTS
if plot_sgrams
    print("Making sgram plots....\n")
    # define timing points for TC+DRUM
    stime_TC = Dates.DateTime(
        Dates.year(TC_DRUM_T[1]),
        Dates.month(TC_DRUM_T[1]),
        Dates.day(TC_DRUM_T[1]),
        Dates.hour(TC_DRUM_T[1])+1
    )
    etime_TC = Dates.DateTime(
        Dates.year(TC_DRUM_T[end]),
        Dates.month(TC_DRUM_T[end]),
        Dates.day(TC_DRUM_T[end]),
        Dates.hour(TC_DRUM_T[end])+1
    )
    time_TC = stime_TC:Dates.Hour(1):etime_TC
    tidx_TC = map(x-> argmin(abs.(TC_DRUM_T .- time_TC[x])), 1:lastindex(time_TC))
    # plot for TC+DRUM
    global hp_sgram_TC = plot(
        TC_DRUM_X,
        TC_DRUM_Y,
        legend = false,
        yflip = true,
        size = (9400,3000),
    )
    annotate!(
        hp_sgram_TC, 
        TC_DRUM_X[tidx_TC], 
        TC_DRUM_Y[tidx_TC] .- 50,
        text.(Dates.format.(time_TC, "m/d-HH"), :black, :left, 20)
    )
    savefig(hp_sgram_TC, string(data_out,"sgram_TC_",suffix,".pdf"))
    # define timing points for PE+DRUM
    stime_PE = Dates.DateTime(
        Dates.year(PE_DRUM_T[1]),
        Dates.month(PE_DRUM_T[1]),
        Dates.day(PE_DRUM_T[1]),
        Dates.hour(PE_DRUM_T[1])+1
    )
    etime_PE = Dates.DateTime(
        Dates.year(PE_DRUM_T[end]),
        Dates.month(PE_DRUM_T[end]),
        Dates.day(PE_DRUM_T[end]),
        Dates.hour(PE_DRUM_T[end])+1
    )
    time_PE = stime_PE:Dates.Hour(1):etime_PE
    tidx_PE = map(x-> argmin(abs.(PE_DRUM_T .- time_PE[x])), 1:lastindex(time_PE))
    # plot for TC+DRUM
    global hp_sgram_PE = plot(
        PE_DRUM_X,
        PE_DRUM_Y,
        legend = false,
        yflip = true,
        size = (9400,3000),
    )
    annotate!(
        hp_sgram_PE, 
        PE_DRUM_X[tidx_PE], 
        PE_DRUM_Y[tidx_PE] .- 50,
        text.(Dates.format.(time_PE, "m/d-HH"), :black, :left, 20)
    )
    savefig(hp_sgram_PE, string(data_out,"sgram_PE_",suffix,".pdf"))
end

## INTERPOLATE EVERYTHING TO COARSEST SAMPLING
print(string("Interpolating everything to coarsest sampling...\n"))
# get target millisecond TC+Q330 timing
global mintmp = maximum([TC_DRUM_T[1], PE_DRUM_T[1], PE_Q330_T[1]])
maxtmp = minimum([TC_DRUM_T[end], PE_DRUM_T[end], PE_Q330_T[end]])
MASTER_dt = maximum([
    TC_Q330_dt, PE_Q330_dt,
    round(mean(Dates.value.(diff(TC_DRUM_T))))/1000,
    round(mean(Dates.value.(diff(PE_DRUM_T))))/1000
])
MASTER_T = mintmp:Dates.Millisecond(convert(Int,round(MASTER_dt*1000))):maxtmp
# make interpolants
function make_interp_tmp(T,D,newT,newdt)
    # locate NaNs relative to newT
    dt = mean(Dates.value.(diff(T)))./1000
    samp_idx = 1:convert(Int,ceil(0.75*newdt/dt)):length(D)
    nan_idx = samp_idx[findall(isnan.(D[samp_idx]))]
    nan_idx = nan_idx[findall(newT[1] .<= T[nan_idx] .<= newT[end])] # get only NaN within MASTER_T
    new_nan_idx = []
    for i = 1:lastindex(nan_idx)
        tmp_idx = argmin(abs.(newT .- T[nan_idx[i]])) # closest newT
        append!(new_nan_idx,tmp_idx) 
    end
    unique!(new_nan_idx) # get unique
    # eliminate nans
    T = T[findall(.!isnan.(D))]
    D = D[findall(.!isnan.(D))]
    if length(T)>length(unique(T)) # if duplicates
        ctmp = countmap(T)
        repeat_idx = findall(values(ctmp).>1)
        duplicate_t = collect(keys(ctmp))[repeat_idx]
        for i = 1:lastindex(duplicate_t)
            rpt_idx = findall(T.==duplicate_t[i])
            D[rpt_idx[1]] = mean(D[rpt_idx])
            deleteat!(T,rpt_idx[2:end])
            deleteat!(D,rpt_idx[2:end])
        end
    end
    # sort
    sidx = sortperm(T)
    T = T[sidx]
    D = D[sidx]
    # get interpolant
    itp = LinearInterpolation(Dates.value.(T .- newT[1]), D)
    # interpolate 
    newD = itp(Dates.value.(newT .- newT[1]))
    # reinsert NaN
    newD[new_nan_idx] .= NaN
    return newD
end
# interpolate
print("0\n")
TC_DRUM_D = make_interp_tmp(TC_DRUM_T,TC_DRUM_D,MASTER_T,MASTER_dt)
print("1\n")
TC_Q330_D = make_interp_tmp(TC_Q330_T,TC_Q330_D,MASTER_T,MASTER_dt)
print("2\n")
PE_DRUM_D = make_interp_tmp(PE_DRUM_T,PE_DRUM_D,MASTER_T,MASTER_dt)
print("3\n")
PE_Q330_D = make_interp_tmp(PE_Q330_T,PE_Q330_D,MASTER_T,MASTER_dt)
print("4\n")

## BANDPASS FILTER
if Band_Filt_All
    print(string("Filtering all data between ",b1_all," and ",b2_all,"Hz\n\n"))
    function tmp_bandpass(T,D,dt) # bandpass function
        # get good indices
        gidx = findall(.!isnan.(D))
        # trim ends if they are NaN
        D = D[gidx[1]:gidx[end]]
        T = T[gidx[1]:gidx[end]]
        # get the gidx again after trimming
        gidx = findall(.!isnan.(D))
        # get NaN
        nan_idx = findall(isnan.(D))
        if !isempty(nan_idx)
        # interpolate if there are NaN
            no_nan_itp = LinearInterpolation(Dates.value.(T[gidx].-T[1]),D[gidx])
            D = no_nan_itp(Dates.value.(T.-T[1]))
        end
        # filter
        tmpfilt = DSP.digitalfilter(
            DSP.Bandpass(b1_all, b2_all; fs = 1/dt),designmethod_all)
        newD = DSP.filt(tmpfilt, D)
        # reintroduce NaN
        if !isempty(nan_idx)
            newD[nan_idx] .= NaN
        end
        return newD, tmpfilt
    end
    # apply filter
    TC_Q330_D, tmpfilt = tmp_bandpass(MASTER_T,TC_Q330_D,MASTER_dt)
    TC_DRUM_D, tmpfilt = tmp_bandpass(MASTER_T,TC_DRUM_D,MASTER_dt)
    PE_Q330_D, tmpfilt = tmp_bandpass(MASTER_T,PE_Q330_D,MASTER_dt)
    PE_DRUM_D, tmpfilt = tmp_bandpass(MASTER_T,PE_DRUM_D,MASTER_dt)
    # plot filter response
    w = range(0,(0.5/TC_Q330_dt),length=250) # get w vector to Nyquist since response will be normalized
    fresp = DSP.freqz(tmpfilt) # frequency response
    presp = DSP.phasez(tmpfilt) # phase response
    hpresp = plot(
        w,
        abs.(fresp),
        title = "Filter Response",
        label = "Freq",
        ylabel = "Freq Response",
        lc = :blue,
    )
    plot!(
        twinx(hpresp),
        w,
        presp,
        label = "Phase",
        ylabel ="Phase Response",
        lc = :red,
        legend = :bottomright,
    )
    plot!(
        hpresp,
        right_margin = 10mm,
    )
    savefig(hpresp,string(data_out,"Bandpass_All_response_",suffix,".pdf"))
end

## CALCULATION OF SMALL WINDOW SCALINGS
print("Calculating Small Window Scalings...\n")
# get small windows
sml_wndw_strt = MASTER_T[1]:sml_wndw_step:MASTER_T[end]-sml_wndw_size
# set up scaling variables (relative to TC+Q330)
TC_Q330_scl = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
TC_Q330_amp = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
TC_DRUM_scl = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
TC_DRUM_amp = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
PE_Q330_scl = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
PE_Q330_amp = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
PE_DRUM_scl = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
PE_DRUM_amp = fill!(Vector{Float64}(undef,length(sml_wndw_strt)),NaN)
# function
function get_scl_and_amp(D,Dref)
    # check for NaN first
    if minimum([sum(.!isnan.(Dref))/length(Dref), sum(.!isnan.(D))/length(D)])>0.9 # nan thresh
        # detrend
        refamp = mean(abs.(filter(!isnan,Dref).-mean(filter(!isnan,Dref))))
        amp = mean(abs.(filter(!isnan,D).-mean(filter(!isnan,D))))
        # compare
        scl =  refamp / amp
    else
        amp = NaN
        scl = NaN
    end
    return  amp, scl
end
# loop over
for i in ProgressBar(1:lastindex(sml_wndw_strt))
    tidx = findall(sml_wndw_strt[i] .<= MASTER_T .<= sml_wndw_strt[i]+sml_wndw_size)
    TC_Q330_amp[i], TC_Q330_scl[i] = get_scl_and_amp(TC_Q330_D[tidx],TC_Q330_D[tidx])
    TC_DRUM_amp[i], TC_DRUM_scl[i] = get_scl_and_amp(TC_DRUM_D[tidx],TC_Q330_D[tidx])
    PE_Q330_amp[i], PE_Q330_scl[i] = get_scl_and_amp(PE_Q330_D[tidx],TC_Q330_D[tidx])
    PE_DRUM_amp[i], PE_DRUM_scl[i] = get_scl_and_amp(PE_DRUM_D[tidx],TC_Q330_D[tidx])
end

## SMALL PLOTS IF DESIRED
if plot_small_segs
    print("Making Small Segment Plots...\n")
    small_seg_strts = MASTER_T[1]:small_seg_step:MASTER_T[end]-small_seg_wndw
    if !isdir(string(data_out,"small_segs_",suffix,"/"))
        mkdir(string(data_out,"small_segs_",suffix,"/"))
    end
    for i in ProgressBar(1:lastindex(small_seg_strts))
        # get indices for sml_wndw_strt and MASTER_T
        idx_MST = findall(small_seg_strts[i] .<= MASTER_T .<= small_seg_strts[i]+small_seg_wndw)
        idx_sml = findall(small_seg_strts[i] .<= sml_wndw_strt.+(0.5*sml_wndw_size) .<= small_seg_strts[i]+small_seg_wndw)
        # plot the traces and scaling values
        hp_time_TC_Q330 = plot(
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            TC_Q330_amp[idx_sml],
            title = "TC+Q330",
            label = "Amplitude",
            ylabel = "Amplitude",
            lc = :blue,
            legend = :topleft,
        )
        plot!(
            hp_time_TC_Q330,
            MASTER_T[idx_MST],
            TC_Q330_D[idx_MST],
            label = "trace",
            lc = :black,
        )
        plot!(
            twinx(hp_time_TC_Q330),
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            TC_Q330_scl[idx_sml],
            label = "Scalar",
            ylabel ="Scalar",
            lc = :red,
            legend = :topright,
        )
        hp_time_TC_DRUM = plot(
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            TC_DRUM_amp[idx_sml],
            title = "TC+DRUM",
            label = "Amplitude",
            ylabel = "Amplitude",
            lc = :blue,
            legend = :topleft,
        )
        plot!(
            hp_time_TC_DRUM,
            MASTER_T[idx_MST],
            TC_DRUM_D[idx_MST],
            label = "trace",
            lc = :black,
        )
        plot!(
            twinx(hp_time_TC_DRUM),
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            TC_DRUM_scl[idx_sml],
            label = "Scalar",
            ylabel ="Scalar",
            lc = :red,
            legend = :topright,
        )
        hp_time_PE_Q330 = plot(
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            PE_Q330_amp[idx_sml],
            title = "PE+Q330",
            label = "Amplitude",
            ylabel = "Amplitude",
            lc = :blue,
            legend = :topleft,
        )
        plot!(
            hp_time_PE_Q330,
            MASTER_T[idx_MST],
            PE_Q330_D[idx_MST],
            label = "trace",
            lc = :black,
        )
        plot!(
            twinx(hp_time_PE_Q330),
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            PE_Q330_scl[idx_sml],
            label = "Scalar",
            ylabel ="Scalar",
            lc = :red,
            legend = :topright,
        )
        hp_time_PE_DRUM = plot(
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            PE_DRUM_amp[idx_sml],
            title = "PE+DRUM",
            label = "Amplitude",
            ylabel = "Amplitude",
            lc = :blue,
            legend = :topleft,
        )
        plot!(
            hp_time_PE_DRUM,
            MASTER_T[idx_MST],
            PE_DRUM_D[idx_MST],
            label = "trace",
            lc = :black,
        )
        plot!(
            twinx(hp_time_PE_DRUM),
            sml_wndw_strt[idx_sml] .+ (0.5*sml_wndw_size),
            PE_DRUM_scl[idx_sml],
            label = "Scalar",
            ylabel ="Scalar",
            lc = :red,
            legend = :topright,
        )
        hp_time_all = plot(
            hp_time_TC_Q330,
            hp_time_TC_DRUM,
            hp_time_PE_Q330,
            hp_time_PE_DRUM,
            layout = grid(4,1),
            size = (1000,1200),
            right_margin = 25mm,
        )
        # save
        savefig(hp_time_all, string(data_out,"small_segs_",suffix,"/",small_seg_strts[i],".pdf"))
    end
end

## PLOT SCALING ON SEISMOGRAMS
if plot_sgrams
    print("Plotting scalings on the seismograms...\n")
    # get Z for each small_seg_strt
    tmp_idx = map(x -> argmin(abs.(Dates.value.(TC_DRUM_T .- (sml_wndw_strt[x] + 0.5*sml_wndw_size)))),
        1:lastindex(sml_wndw_strt)
    )
    t_tmp = TC_DRUM_T[tmp_idx]
    x_tmp = TC_DRUM_X[tmp_idx]
    z_tmp = TC_DRUM_Z[tmp_idx]
    # add the scl (demedian and normalize first)
    DC_tmp = median(filter(!isnan,TC_DRUM_scl))
    scl_tmp = TC_DRUM_scl .- DC_tmp
    Y_rng_tmp = maximum(filter(!isnan,TC_DRUM_Y))-minimum(filter(!isnan,TC_DRUM_Y))
    amp_tmp = maximum(abs.(filter(!isnan,TC_DRUM_D)))/maximum(abs.(filter(!isnan,scl_tmp)))
    scl_tmp = scl_tmp .* -amp_tmp
    y_tmp = z_tmp .+ scl_tmp
    # insert NaN 
    nan_idx_0 = findall(isnan,TC_DRUM_X)
    nan_idx = map(x-> argmin(abs.(Dates.value.(t_tmp .- TC_DRUM_T[nan_idx_0[x]]))),
        1:lastindex(nan_idx_0)
    )
    y_tmp[nan_idx] .= NaN
    x_tmp[nan_idx] .= NaN
    # plot for TC
    plot!(
        hp_sgram_TC,
        x_tmp,
        y_tmp,
        lc = :red,
    )
    # save for TC
    savefig(hp_sgram_TC,string(data_out,"sgram_scl_TC_",suffix,".pdf"))
    # get Z for each small_seg_strt
    tmp_idx = map(x -> argmin(abs.(Dates.value.(PE_DRUM_T .- (sml_wndw_strt[x] + 0.5*sml_wndw_size)))),
        1:lastindex(sml_wndw_strt)
    )
    t_tmp = PE_DRUM_T[tmp_idx]
    x_tmp = PE_DRUM_X[tmp_idx]
    z_tmp = PE_DRUM_Z[tmp_idx]
    # add the scl (demedian and normalize first)
    DC_tmp = median(filter(!isnan,PE_DRUM_scl))
    scl_tmp = PE_DRUM_scl .- DC_tmp
    Y_rng_tmp = maximum(filter(!isnan,PE_DRUM_Y))-minimum(filter(!isnan,PE_DRUM_Y))
    amp_tmp = maximum(abs.(filter(!isnan,PE_DRUM_D)))/maximum(abs.(filter(!isnan,scl_tmp)))
    scl_tmp = scl_tmp .* -amp_tmp
    y_tmp = z_tmp .+ scl_tmp
    # insert NaN 
    nan_idx_0 = findall(isnan,PE_DRUM_X)
    nan_idx = map(x-> argmin(abs.(Dates.value.(t_tmp .- PE_DRUM_T[nan_idx_0[x]]))),
        1:lastindex(nan_idx_0)
    )
    y_tmp[nan_idx] .= NaN
    x_tmp[nan_idx] .= NaN
    # plot for PE
    plot!(
        hp_sgram_PE,
        x_tmp,
        y_tmp,
        lc = :red,
    )
    # save for PE
    savefig(hp_sgram_PE,string(data_out,"sgram_scl_PE_",suffix,".pdf"))
end

print("Making Summary Plots....\n")
## PLOT HISTOGRAMS
hp_hist_TC_DRUM_amp = histogram(
    TC_DRUM_amp,
    title = "TC+DRUM amp",
    legend = false,
)
hp_hist_TC_DRUM_scl = histogram(
    TC_DRUM_scl,
    title = "TC+DRUM scl",
    legend = false,
)
hp_hist_PE_Q330_amp = histogram(
    PE_Q330_amp,
    title = "PE+Q330 amp",
    legend = false,
)
hp_hist_PE_Q330_scl = histogram(
    PE_Q330_scl,
    title = "PE+Q330 scl",
    legend = false,
)
hp_hist_PE_DRUM_amp = histogram(
    PE_DRUM_amp,
    title = "PE+DRUM amp",
    legend = false,
)
hp_hist_PE_DRUM_scl = histogram(
    PE_DRUM_scl,
    title = "PE+DRUM scl",
    legend = false,
)
hp_hist_all = plot(
    hp_hist_TC_DRUM_amp,
    hp_hist_PE_Q330_amp,
    hp_hist_PE_DRUM_amp,
    hp_hist_TC_DRUM_scl,
    hp_hist_PE_Q330_scl,
    hp_hist_PE_DRUM_scl,
    layout = grid(2,3),
    size = (1200,800),
)
savefig(hp_hist_all, string(data_out,"histogram_",suffix,".pdf"))

## PLOT WITH TIME
hp_time_TC_Q330 = plot(
    sml_wndw_strt,
    TC_Q330_amp,
    title = "TC+Q330",
    label = "Amplitude",
    ylabel = "Amplitude",
    lc = :blue,
    legend = :topleft,
)
plot!(
    twinx(hp_time_TC_Q330),
    sml_wndw_strt,
    TC_Q330_scl,
    label = "Scalar",
    ylabel ="Scalar",
    lc = :red,
    legend = :topright,
)
hp_time_TC_DRUM = plot(
    sml_wndw_strt,
    TC_DRUM_amp,
    title = "TC+DRUM",
    label = "Amplitude",
    ylabel = "Amplitude",
    lc = :blue,
    legend = :topleft,
)
plot!(
    twinx(hp_time_TC_DRUM),
    sml_wndw_strt,
    TC_DRUM_scl,
    label = "Scalar",
    ylabel ="Scalar",
    lc = :red,
    legend = :topright,
)
hp_time_PE_Q330 = plot(
    sml_wndw_strt,
    PE_Q330_amp,
    title = "PE+Q330",
    label = "Amplitude",
    ylabel = "Amplitude",
    lc = :blue,
    legend = :topleft,
)
plot!(
    twinx(hp_time_PE_Q330),
    sml_wndw_strt,
    PE_Q330_scl,
    label = "Scalar",
    ylabel ="Scalar",
    lc = :red,
    legend = :topright,
)
hp_time_PE_DRUM = plot(
    sml_wndw_strt,
    PE_DRUM_amp,
    title = "PE+DRUM",
    label = "Amplitude",
    ylabel = "Amplitude",
    lc = :blue,
    legend = :topleft,
)
plot!(
    twinx(hp_time_PE_DRUM),
    sml_wndw_strt,
    PE_DRUM_scl,
    label = "Scalar",
    ylabel ="Scalar",
    lc = :red,
    legend = :topright,
)
hp_time_all = plot(
    hp_time_TC_Q330,
    hp_time_TC_DRUM,
    hp_time_PE_Q330,
    hp_time_PE_DRUM,
    layout = grid(4,1),
    size = (1000,1200),
    right_margin = 25mm,
)
savefig(hp_time_all, string(data_out,"timeseries_",suffix,".pdf"))

## PLOT SCATTERS
hp_TC_DRUM = scatter(
    TC_DRUM_amp,
    TC_DRUM_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="TC+DRUM",
    legend=false,
)
hp_PE_Q330 = scatter(
    PE_Q330_amp,
    PE_Q330_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="PE+Q330",
    legend=false,
)
hp_PE_DRUM = scatter(
    PE_DRUM_amp,
    PE_DRUM_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="PE+DRUM",
    legend=false,
)
cidx = findall(
    sort(TC_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= TC_DRUM_amp .<= 
    sort(TC_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_TC_DRUM_cull = scatter(
    TC_DRUM_amp[cidx],
    TC_DRUM_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("TC+DRUM (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
)
cidx = findall(
    sort(PE_Q330_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= PE_Q330_amp .<= 
    sort(PE_Q330_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_PE_Q330_cull = scatter(
    PE_Q330_amp[cidx],
    PE_Q330_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("PE+Q330 (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
)
cidx = findall(
    sort(PE_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= PE_DRUM_amp .<= 
    sort(PE_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_PE_DRUM_cull = scatter(
    PE_DRUM_amp[cidx],
    PE_DRUM_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("PE+DRUM (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
)
hp_all = plot(
    hp_TC_DRUM,
    hp_PE_Q330,
    hp_PE_DRUM,
    hp_TC_DRUM_cull,
    hp_PE_Q330_cull,
    hp_PE_DRUM_cull,
    layout=grid(2,3),
    size=(1400,900)
)
savefig(hp_all,string(data_out,"scaling_factors_",suffix,".pdf"))

## PLOT SCATTERS IN LOG LOG
hp_TC_DRUM = scatter(
    TC_DRUM_amp,
    TC_DRUM_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="TC+DRUM",
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
hp_PE_Q330 = scatter(
    PE_Q330_amp,
    PE_Q330_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="PE+Q330",
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
hp_PE_DRUM = scatter(
    PE_DRUM_amp,
    PE_DRUM_scl,
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title="PE+DRUM",
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
cidx = findall(
    sort(TC_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= TC_DRUM_amp .<= 
    sort(TC_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_TC_DRUM_cull = scatter(
    TC_DRUM_amp[cidx],
    TC_DRUM_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("TC+DRUM (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
cidx = findall(
    sort(PE_Q330_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= PE_Q330_amp .<= 
    sort(PE_Q330_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_PE_Q330_cull = scatter(
    PE_Q330_amp[cidx],
    PE_Q330_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("PE+Q330 (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
cidx = findall(
    sort(PE_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*lower_percentile))+1]
    .<= PE_DRUM_amp .<= 
    sort(PE_DRUM_amp)[convert(Int,round(length(sml_wndw_strt)*upper_percentile))])
hp_PE_DRUM_cull = scatter(
    PE_DRUM_amp[cidx],
    PE_DRUM_scl[cidx],
    ma = 0.3,
    xlabel="amplitude",
    ylabel="scaling",
    title=string("PE+DRUM (",lower_percentile*100,"-",upper_percentile*100,"%)"),
    legend=false,
    yaxis=:log,
    xaxis=:log,
)
hp_all = plot(
    hp_TC_DRUM,
    hp_PE_Q330,
    hp_PE_DRUM,
    hp_TC_DRUM_cull,
    hp_PE_Q330_cull,
    hp_PE_DRUM_cull,
    layout=grid(2,3),
    size=(1400,900)
)
savefig(hp_all,string(data_out,"log_scaling_factors_",suffix,".pdf"))

print("Done!")