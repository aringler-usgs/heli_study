# compare_xy_sac.jl
#= 
This script will read in DigitSeis xy data as a csv and then compute the 
time differences between it and relevant SAC files by running xcorr. The
code will output a map of timing uncertainties / offsets as well as areas
where data is missing.

Created on: 11/02/2022
Written by: Thomas Lee

Last Modified: 11/14/2022
    - added option to filter sac
    - added first dericative of tdiff plot
    - added automated setup text write out
    - flipped xy and sac in calculation so that Tdiff reflects the correction xy needs
Last Modified: 11/22/2022
    - added exit case for NaN in drift cross-corr correction
Last Modified: 2/2/2023
    - added writeout of original x-y positions with associated times
Last Modified: 2/21/2023
    - changed making y offset the same between traces to subtracting a difference of the
      average y positions from one trace to the next
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
using Distributions
import DSP
using Measures

## SETUP
# csv file
csv_in = "/Users/thomaslee/Research/Microform/ASL_NewRecords/Digitized/xy_new/DigitSeis_20220714_TCZ_ASL.trace.txt"
# sac or dir of sac
sac_in = "/Users/thomaslee/Research/Microform/ASL_NewRecords/ASL_SAC/00_TCZ/"
# output data and plots
# cDataOut = "/Users/thomaslee/Research/Microform/ASL_NewRecords/XY_timeseries/"
cDataOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220714/"
cPlotOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220714/"
cDataOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_test/"
cPlotOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_test/"
# cDataOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220713_test/"
# cPlotOut = "/Users/thomaslee/Desktop/XY_SAC_RESULTS_TCZ_20220713_test/"
# parameters for detrending
trim_trace = 4
dtrnd_wndw = 2500 # n points in detrending window
dtrnd_step = 100 # n points between detrending windows
# estimated start and end times
# est_stime = Dates.DateTime(2022,7,13,21,23) # estimated start time of xy
# est_etime = Dates.DateTime(2022,7,14,14,25) # estimated end time of xy
est_stime = Dates.DateTime(2022,7,14,14,25) # estimated start time of xy
est_etime = Dates.DateTime(2022,7,20,17,5) # estimated end time of xy
est_error = Dates.Second(600) # error on estimate (overestimate since this determines width of search)
est_error_step = Dates.Second(1) # accuracy to which to nail get start and end times
est_strt_end_wndw = Dates.Second(300) # size of window to look for
stime_off = Dates.Second(300) # time after estimated start time to throw out
etime_off = Dates.Second(10) # time before estimated end time to throw out
fit_thresh = 0.9 # throw warning if below this value
# cross corr parameters for drift]]
ccstep = Dates.Second(10)
ccwndw = Dates.Second(180) # must be even!!
ccrang = Dates.Millisecond.(-20000:250:20000)
weightCCdrift = true
# plot parameters
plot_cc_diag = true
plot_small_segs = true
small_seg_wndw = Dates.Minute(15)
# filter parameters
useFilt = false # filter the sac to emulate an analog filter
b1 = 0.01 # lowband
b2 = 1.34 # highband
designmethod = DSP.Butterworth(2)

## SETUP DIRECTORIES
# plots
if !isdir(cDataOut)
    mkdir(cDataOut)
end
# data output
if !isdir(cPlotOut)
    mkdir(cPlotOut)
end

## CSV READ IN
# get xy_name
xy_name = csv_in[(findlast("/",csv_in)[1]+1):end] # get just filename
xy_name = xy_name[1:findfirst(".",xy_name)[1]-1] # get rid of ".trace.txt"
# open file
f = open(csv_in)
# read file
ln = readlines(f)
# close file
close(f)
# get number of lines and length and width
tmpln = split(ln[3],' ') # delimit with spaces
ntraces = parse(Int,tmpln[1])
im_hght = parse(Int,tmpln[2])
im_lgth = parse(Int,tmpln[3])
# intialize trace variables
trace_x = Vector{Vector{Float32}}(undef,ntraces)
trace_y = Vector{Vector{Float32}}(undef,ntraces)
# loop over lines
global i_line = 4 # skip the three headers
while i_line<=length(ln)
    # check if info for new trace is starting
    if isempty(ln[i_line]) && i_line<length(ln) # trace separator
        global i_line += 1 # move i_line
        lntmp = split(ln[i_line],' ') # delimit with spaces
        global i_trace = parse(Int,lntmp[1]) # get trace number
        trace_x[i_trace] = parse(Int,lntmp[2]):parse(Int,lntmp[3]) # get x-values
        global n_pts = parse(Int,lntmp[4]) # get npoints in trace
        # check that values align
        if n_pts != length(trace_x[i_trace])
            error(string("On line ",i_line," expected # of x values does not correspond to start and end values!\n"))    
        end
        # get next n_pts to read in trace
        trace_y[i_trace] = map(x->parse(Float32,ln[i_line+x]),1:n_pts)
        # update i_line to next value
        i_line = i_line + n_pts + 1
    elseif isempty(ln[i_line]) && i_line==length(ln) # file end
        if i_trace==ntraces # check trace numbers
            i_line += 1 # update i_line to exit loop
        else
            error("Number of traces read (",i_trace,") does not match the number expected (",ntraces,")!\n")
        end
    else
        error(string("Expected a trace seperator (empty line) at line ",i_line,", but got: ",ln[i_line],"!\n"))
    end
end
# report
print(string("Read in ",ntraces," traces from ",xy_name,".trace.txt\n"))
# perform stitiching 
for i = 1:lastindex(trace_x)
    if i==1 # initialize globals
        global xy_X = trace_x[i].-(trace_x[i][1]-1) # set first x to 0
        global xy_Y = copy(trace_y[i]) 
    else # append with y-matching and x-adding
        # get data
        tmp_X = trace_x[i][trim_trace:end-trim_trace]
        tmp_Y = trace_y[i][trim_trace:end-trim_trace]
        # take out DC offset
        # tmp_Y = tmp_Y .- (tmp_Y[1]-xy_Y[end])
        # tmp_Y = tmp_Y .- (median(filter(!isnan,trace_y[i]))-median(filter(!isnan,trace_y[i-1]))) # average subtraction (doesn't work since its left end to right end)
        tmp_Y = tmp_Y .- (tmp_Y[1]-xy_Y[end]-(((xy_Y[end]-xy_Y[end-1]) + (tmp_Y[2]-tmp_Y[1]))/2)) # average of rate of change on either side
        # make X continuous
        tmp_X = tmp_X .- (tmp_X[1]-xy_X[end]-1)
        # append
        append!(xy_X, tmp_X)
        append!(xy_Y, tmp_Y)
    end
end
print(string("Stitched traces together to get ",length(xy_Y)," samples\n"))
# take out polarity flip (down in digitseis is up on seismogram)
xy_Y = -1 .* xy_Y
xy_Y0 = copy(xy_Y) # save original before any detrending
# take out linear trend (using start and ends of traces)
strt_end_idx = vcat((1:dtrnd_wndw),(length(xy_Y)-dtrnd_wndw:length(xy_Y)))
x = xy_X[strt_end_idx[findall(!isnan,xy_Y[strt_end_idx])]]
y = xy_Y[strt_end_idx[findall(!isnan,xy_Y[strt_end_idx])]]
a = ((length(x)*sum(x.*y)-(sum(x)*sum(y))) /  # solve for slope
    ((length(x)*sum(x.^2))-(sum(x)^2)))
b = (1/(length(x)))*(sum(y)-a*sum(x)) # solve for y-intercept
trend_y_linear = a*xy_X .+ b # get trend
# subtract trend from y0
xy_Y = xy_Y .- trend_y_linear
xy_Y1 = copy(xy_Y) # save linear version
# run non-parametric detrending
dtrnd_strt = convert.(Int,xy_X[1]:dtrnd_step:(xy_X[end]-dtrnd_wndw)) # get starts (for computing values)
dtrnd_ctr = convert.(Int,round.(dtrnd_strt .+ dtrnd_wndw/2)) # get centers (for interpolating)
# get mean of windows
dtrnd_y = map(x->
    mean(filter(!isnan,xy_Y[dtrnd_strt[x]:(dtrnd_strt[x]+dtrnd_wndw)])),
    1:lastindex(dtrnd_strt),
) 
# interpolate non-parametric line and then subtract
itp_y = LinearInterpolation(dtrnd_ctr, dtrnd_y)
trend_y = itp_y(dtrnd_ctr[1]:dtrnd_ctr[end]) # may need to remove NaN and then put back in later
# taper trend and pad with zeros
trend_y_tapered = vcat( 
    zeros(dtrnd_ctr[1]-1),
    trend_y.*vcat( # cosine taper with flat top
        (cos.(range(pi,2*pi,length=dtrnd_ctr[1])).+1)./2, # from -1 to 1
        ones(length(trend_y)-2*dtrnd_ctr[1]),
        (cos.(range(0,pi,length=dtrnd_ctr[1])).+1)./2
    ), # from 1 to -1 
    zeros(length(xy_X)-dtrnd_ctr[end]),
) 
# get detrended version
xy_Y = xy_Y .- trend_y_tapered

## PLOT XY PRODUCTS
# plot seismogram style
hp_sgram = plot(
    trace_x[1],
    trace_y[1],
    legend = false,
    title = xy_name,
    size = (im_lgth,im_hght),
    yflip = true,
)
for i = 2:lastindex(trace_x)
    plot!(
        hp_sgram,
        trace_x[i],
        trace_y[i],
    )
end
savefig(hp_sgram,string(cPlotOut,xy_name,"_sgram.pdf"))
# plot detrended and stitched
hp_orig = plot(
    xy_X,
    xy_Y0,
    lc = :blue,
    title = xy_name,
    label = "original",
)
plot!(
    hp_orig,
    xy_X,
    trend_y_linear,
    lc = :red,
    label = "linear trend",
)
hp_ord1 = plot(
    xy_X,
    xy_Y1,
    lc = :blue,
    label = "linear trend removed",
)
plot!(
    hp_ord1,
    xy_X,
    trend_y_tapered,
    lc = :red,
    label = "rolling avg. trend",
)
hp_npar = plot(
    xy_X,
    xy_Y,
    lc = :blue,
    label = "non-parametric trend removed",
)
hp_trend = plot(
    hp_orig,
    hp_ord1,
    hp_npar,
    layout = grid(3,1),
    size = (1000,1000),
)
savefig(hp_trend,string(cPlotOut,xy_name,"_detrend.pdf"))
print(string("Saved xy_plots to ",cPlotOut,"\n"))

## SAVE XY DETRENDED AND STITCHED 
f = open(string(cDataOut,xy_name,".txt"),"w")
print(f,string("stitched_x,stitched_y,detrend_y,\n"))
for i = 1:lastindex(xy_X)
    print(f,string(xy_X[i],",",xy_Y0[i],",",xy_Y[i],",\n"))
end
close(f)
print(string("Saved detrended and stitiched xy_data to ",cDataOut,xy_name,".txt\n"))
print("\n")

## SAC READ IN
print(string("Reading sac from ",sac_in,"\n"))
# check if is dir or sac
Stmp = lf.combsac(sac_in)
names = string(Stmp.sta.sta,".",Stmp.sta.cha)
sac_dt = round(Stmp.delta,digits=2)
sac_T = lf.gettime(Stmp)
sac_D = trace(Stmp)
print(string("Read in ",length(sac_D)," samples\n"))
print("\n")

## FILTER IF DESIRED
if useFilt
    print(string("Filtering the data between ",b1," and ",b2,"Hz\n\n"))
    # get good indices
    local gidx = findall(.!isnan.(sac_D))
    # trim ends if they are NaN
    sac_D = sac_D[gidx[1]:gidx[end]]
    sac_T = sac_T[gidx[1]:gidx[end]]
    # get tthe gidx again after trimming
    gidx = findall(.!isnan.(sac_D))
    # get NaN
    nan_idx = findall(isnan.(sac_D))
    # interpolate
    no_nan_itp = LinearInterpolation(
        Dates.value.(sac_T[gidx].-sac_T[1]),
        sac_D[gidx]
    )
    sac_D = no_nan_itp(Dates.value.(sac_T.-sac_T[1]))
    # filter
    tmpfilt = DSP.digitalfilter(
        DSP.Bandpass(b1, b2; fs = 1/sac_dt),
        designmethod,
    )
    sac_D = DSP.filt(
        tmpfilt, 
        sac_D,
    )
    # reintroduce NaN
    sac_D[nan_idx] .= NaN
    # plot filter response
    w = range(0,(0.5/sac_dt),length=250) # get w vector to Nyquist since response will be normalized
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
    savefig(hpresp,string(cPlotOut,xy_name,"_response.pdf"))
end

## CROSS CORRELATE START AND END
print("Searching for best fitting start and end times\n")
# estimate sample rate from length of xy and estimated times
est_samp_rate = length(xy_X)/(Dates.value(est_etime-est_stime)/1000)
# get the start and end snippets from the xy_data
if typeof.(est_strt_end_wndw).==Dates.Second
    tmp_len = convert(Int,round(Dates.value(est_strt_end_wndw)*est_samp_rate))
    tmp_soff = convert(Int,round(Dates.value(stime_off)*est_samp_rate))
    tmp_eoff = convert(Int,round(Dates.value(etime_off)*est_samp_rate))
    global xy_strt_Y = xy_Y[1+tmp_soff:tmp_len+tmp_soff]
    global xy_strt_X = xy_X[1+tmp_soff:tmp_len+tmp_soff]
    global xy_end_Y = xy_Y[end-tmp_len-tmp_eoff:end-tmp_eoff]
    global xy_end_X = xy_X[end-tmp_len-tmp_eoff:end-tmp_eoff]
else
    error(("!!! est_strt_end_wndw, stime_off, and etime_off must be ::Dates.Second !!!\n"))
end
# get windowing of sac within estimated error and of size error
wndws_strt = (est_stime-est_error+stime_off):est_error_step:(est_stime+est_error+stime_off) # start points of start windows
wndws_end = (est_etime-est_error-etime_off):est_error_step:(est_etime+est_error-etime_off) # end points of end windows
# loop over windows and log cc
global cc_coef_start = Vector{Float32}(undef,length(wndws_strt))
global cc_coef_end = Vector{Float32}(undef,length(wndws_end))
global xy_cc_strt = []
global sac_cc_strt = []
global xy_cc_end = []
global sac_cc_end = []
for i in ProgressBar(1:lastindex(wndws_strt))
    # get sac for that time
    strt_idx = findall(wndws_strt[i] .<= sac_T .<= (wndws_strt[i]+est_strt_end_wndw))
    end_idx = findall((wndws_end[i]-est_strt_end_wndw) .<= sac_T .<= wndws_end[i])
    # interpolate data with more samples to the same sampling
    if (1/sac_dt)>=est_samp_rate # interpolate the sac
        # do start
        old_T_strt = Dates.value.(sac_T[strt_idx].-sac_T[strt_idx[1]]) # the integer valued times
        old_T_end = Dates.value.(sac_T[end_idx].-sac_T[end_idx[1]])
        new_T_strt = range(old_T_strt[1],old_T_strt[end],length=length(xy_strt_Y))
        new_T_end = range(old_T_end[1],old_T_end[end],length=length(xy_end_Y))
        itp_sac_strt = LinearInterpolation(old_T_strt, sac_D[strt_idx])
        itp_sac_end = LinearInterpolation(old_T_end, sac_D[end_idx])
        global sac_strt = itp_sac_strt(new_T_strt)
        global sac_end = itp_sac_end(new_T_end)
        if isempty(xy_cc_strt)
            global xy_end = xy_end_Y
            global xy_strt = xy_strt_Y # get the uninterpolate data
        end
        # do end
    else # interpolate the xy
        if isempty(xy_cc_strt)
            new_T_strt = range(xy_strt_X[1],xy_strt_X[end],length=length(strt_idx))
            new_T_end = range(xy_end_X[1],xy_end_X[end],length=length(end_idx))
            itp_xy_strt = LinearInterpolation(xy_strt_X, xy_strt_Y)
            itp_xy_end = LinearInterpolation(xy_end_X, xy_end_Y)
            global xy_strt = itp_xy_strt(new_T_strt)
            global xy_end = itp_xy_end(new_T_end)
        end
        global sac_strt = sac_D[strt_idx] # get the uninterpolated data
        global sac_end = sac_D[end_idx]
    end
    # demean and normalize
    sac_strt = sac_strt .- mean(filter(!isnan,sac_strt))
    sac_strt = sac_strt ./ median(abs.(filter(!isnan,sac_strt)))
    push!(sac_cc_strt,sac_strt)
    sac_end = sac_end .- mean(filter(!isnan,sac_end))
    sac_end = sac_end ./ median(abs.(filter(!isnan,sac_end)))
    push!(sac_cc_end,sac_end)
    if isempty(xy_cc_strt) # add the single set of xy data
        xy_strt = xy_strt .- mean(filter(!isnan,xy_strt))
        xy_strt = xy_strt ./ median(abs.(filter(!isnan,xy_strt)))
        append!(xy_cc_strt,xy_strt)
        xy_end = xy_end .- mean(filter(!isnan,xy_end))
        xy_end = xy_end ./ median(abs.(filter(!isnan,xy_end)))
        append!(xy_cc_end,xy_end)
    end
    # get rid of nan
    gidx_strt = intersect(findall(.!isnan.(sac_strt)),findall(.!isnan.(xy_strt)))
    sac_s = sac_strt[gidx_strt]
    xy_s = xy_strt[gidx_strt]
    gidx_end = intersect(findall(.!isnan.(sac_end)),findall(.!isnan.(xy_end)))
    sac_e = sac_end[gidx_end]
    xy_e = xy_end[gidx_end]
    # get sample cross correlation (normalized)
    cc_coef_start[i] = sum((sac_s .- mean(sac_s)).*(xy_s .- mean(xy_s)))/
        sqrt(sum((sac_s .- mean(sac_s)).^2)*sum((xy_s .- mean(xy_s)).^2))
    cc_coef_end[i] = sum((sac_e .- mean(sac_e)).*(xy_e .- mean(xy_e)))/
        sqrt(sum((sac_e .- mean(sac_e)).^2)*sum((xy_e .- mean(xy_e)).^2))
end
# get best fitting cc
strt_i = argmax(cc_coef_start)
end_i = argmax(cc_coef_end)
xy_t_strt = wndws_strt[strt_i]-stime_off
xy_t_end = wndws_end[end_i]+etime_off
print(string("  Estimated start time is ",Dates.format(est_stime, "yyyy-mm-dd HH:MM:SS"),"\n"))
print(string("  Best fitting start time is ",Dates.format(xy_t_strt, "yyyy-mm-dd HH:MM:SS")," with CCC of ",cc_coef_start[strt_i],"\n"))
if cc_coef_start[strt_i] .< fit_thresh
    print(string("    !!! WARNING: CC value for start (",cc_coef_start[strt_i],") is below threshold (",fit_thresh,") !!!\n"))
end
print(string("  Estimated end time is ",Dates.format(est_etime, "yyyy-mm-dd HH:MM:SS"),"\n"))
print(string("  Best fitting end time is ",Dates.format(xy_t_end, "yyyy-mm-dd HH:MM:SS")," with CCC of ",cc_coef_end[end_i],"\n"))
if cc_coef_end[end_i] .< fit_thresh
    print(string("    !!! WARNING: CC value for end (",cc_coef_end[end_i],") is below threshold (",fit_thresh,") !!!\n"))
end
# plot results and best fit for the start and end
hp_cc_strt = plot( # plot start ccs
    wndws_strt.-stime_off,
    cc_coef_start,
    legend=false,
    title="Start CCs"
)
hp_wndw_strt = plot( # plot best fitting window with 
    [xy_cc_strt, sac_cc_strt[strt_i]],
    label = ["xy" "sac"],
    xlabel = "samples",
    title = string("Best Start Fit ",Dates.format(wndws_strt[strt_i],"yyyy-mm-dd HH:MM:SS")),
) 
hp_start = plot( # aggregate
    hp_cc_strt,
    hp_wndw_strt,
    layout = grid(2,1),
    size = (1500,1500),
)
hp_cc_end = plot( # plot the end ccs
    wndws_end.+etime_off,
    cc_coef_end,
    legend=false,
    title="End CCs"
)
hp_wndw_end = plot( # ending window
    [xy_cc_end, sac_cc_end[end_i]],
    label = ["xy" "sac"],
    xlabel = "samples",
    title = string("Best End Fit ",Dates.format(wndws_end[end_i],"yyyy-mm-dd HH:MM:SS")),
)
hp_end = plot( # aggregate
    hp_cc_end,
    hp_wndw_end,
    layout = grid(2,1),
    size = (1500,1500),
)
# save figs
savefig(hp_start, string(cPlotOut,xy_name,"_start.pdf"))
savefig(hp_end, string(cPlotOut,xy_name,"_end.pdf"))
print(string("Saved figures to ",cPlotOut,"\n\n"))

## LINEARLY INTERPOLATE TIME
xy_lgth = Dates.value(xy_t_end - xy_t_strt) # length in millseconds
xy_T = Dates.Millisecond.(round.(range(0,xy_lgth,length=length(xy_X)))) # evenly space
xy_T = xy_t_strt .+ xy_T # set going forward from start
linear_samp_rate = length(xy_X)/(xy_lgth/1000)

## GET VERSION OF HIGHER SAMPLED DATA THAT IS INTERPOLATED
if linear_samp_rate > (1/sac_dt) # interpolate xy
    tmpidx = findall(xy_t_strt .< sac_T .< xy_t_end)
    global t_int = sac_T[tmpidx]
    global sac_int = sac_D[tmpidx]
    itp = LinearInterpolation(Dates.value.(xy_T.-xy_t_strt),xy_Y)
    global xy_int = itp(Dates.value.(t_int.-xy_t_strt))
else # interpolate sac
    tmpidx = findall((xy_t_strt - Dates.Minute(1)) .< sac_T .< (xy_t_end + Dates.Minute(1)))
    global t_int = xy_T
    global xy_int = xy_Y
    itp = LinearInterpolation(Dates.value.(sac_T[tmpidx].-xy_t_strt),sac_D[tmpidx])
    global sac_int = itp(Dates.value.(t_int.-xy_t_strt))
end

## GET DRIFT IN BETWEEN
# print status
print("Computing drifts from sac to xy\n")
# set up the windows in xy
cc_wind_strt = xy_t_strt:ccstep:(xy_t_end-ccwndw)
cc_wind_ctr = cc_wind_strt .+ (ccwndw/2)
# intialize the vectors of time differences and CC vals
Tdiff = Vector{Dates.Millisecond}(undef,length(cc_wind_strt))
CC = Vector{Float32}(undef,length(cc_wind_strt))
# loop over cc windows
for i in ProgressBar(1:lastindex(cc_wind_strt))
    # get the window for xy
    global xy_tmp = xy_int[findall(cc_wind_strt[i].<= t_int .<= cc_wind_strt[i]+ccwndw)]
    # check if all NaN
    if sum(isnan.(xy_tmp))<length(xy_tmp)
        # demean and normalize
        xy_tmp = xy_tmp .- mean(filter(!isnan,xy_tmp))
        xy_tmp = xy_tmp ./ median(abs.(filter(!isnan,xy_tmp)))
        # set up range of times to try for the time correction needed for xy
        tmp_stimes = cc_wind_strt[i] .+ ccrang
        # set up CC val variable
        tmp_ccs = Vector{Float32}(undef,length(tmp_stimes))
        if plot_cc_diag # save values if diagnostic
            global sac_tmp_segs = []
        end
        # loop over the new xy windows
        for j = 1:lastindex(tmp_stimes)
            # grab sac data in window
            sac_tmp = sac_int[findall(tmp_stimes[j] .<= t_int .<= tmp_stimes[j]+ccwndw)]
            # demean and normalize xy
            sac_tmp = sac_tmp .- mean(filter(!isnan,sac_tmp))
            sac_tmp = sac_tmp ./ median(abs.(filter(!isnan,sac_tmp)))
            # interpolate to number of samples in shorter trace
            if length(xy_tmp)<length(sac_tmp) # interpolate sac
                tmp_itp = LinearInterpolation(
                    range(1,length(xy_tmp),length=length(sac_tmp)),
                    sac_tmp
                )
                global sac_tmp = tmp_itp(1:length(xy_tmp))
                global Btmp = sac_tmp
                global Atmp = xy_tmp
            elseif length(xy_tmp)>length(sac_tmp) # interpolate xy
                tmp_itp = LinearInterpolation(
                    range(1,length(sac_tmp),length=length(xy_tmp)),
                    xy_tmp
                )
                global Atmp = tmp_itp(1:length(sac_tmp))
                global Btmp = sac_tmp
            else
                global Atmp = xy_tmp
                global Btmp = sac_tmp
            end
            # save value if desired
            if plot_cc_diag
                push!(sac_tmp_segs, Btmp)
            end
            # get non-nan values
            gidx = intersect(findall(.!isnan.(Atmp)),findall(.!isnan.(Btmp)))
            Atmp = Atmp[gidx]
            Btmp = Btmp[gidx]
            # calculate CC
            tmp_ccs[j] = sum((Atmp .- mean(Atmp)).*(Btmp.- mean(Btmp)))/
                sqrt(sum((Atmp .- mean(Atmp)).^2)*sum((Btmp .- mean(Btmp)).^2))
        end
        # calculate highest CC
        if weightCCdrift # use weighting function based on CC[i] and Tdiff[i-1]
            if i>1
                # get normal centered on last Tdiff
                global wght_tmp = pdf.(Normal(
                        Dates.value(Tdiff[i-1]),
                        (Dates.value(ccrang[end])-Dates.value(ccrang[1]))/6,
                    ),
                    Dates.value.(ccrang)
                )
            else # center on zero to start
                global wght_tmp = pdf.(Normal(
                        0,
                        (Dates.value(ccrang[end])-Dates.value(ccrang[1]))/6,
                    ),
                    Dates.value.(ccrang)
                )
            end
            # set yrange of pdf to [0,1]
            wght_tmp = wght_tmp .- minimum(wght_tmp)
            wght_tmp = wght_tmp ./ maximum(wght_tmp)
            # set to vary in height based on CC
            wght_tmp = (wght_tmp .* (1-maximum(tmp_ccs))) .+ maximum(tmp_ccs)
            #  now get the highest CC with this weighting
            global i_tmp = argmax(tmp_ccs .* wght_tmp)
        else # simply get the highest cc value and relate to Tdiff   
            global i_tmp = argmax(tmp_ccs)
        end
        # get CC and Tdiff
        CC[i] = tmp_ccs[i_tmp]
        Tdiff[i] = ccrang[i_tmp]
        # plot if we're doing that
        if plot_cc_diag
            if !isdir(string(cPlotOut,xy_name,"_cc_diag/"))
                mkdir(string(cPlotOut,xy_name,"_cc_diag/"))
            end
            hp_tmp_cc = plot(
                Dates.value.(ccrang)./1000,
                tmp_ccs,
                label = "CCs",
                title = string("CCs for ",cc_wind_ctr[i]," best offset corrects xy by ",Dates.value(Tdiff[i])/1000," seconds"),
            ) 
            scatter!(
                hp_tmp_cc,
                [Dates.value(Tdiff[i])/1000],
                [CC[i]],
                mc = :black,
                label = "Max CC",
                shape = :star5,
                ms = 15*CC[i],
            )
            if weightCCdrift # add plots of weighting and weighted CCs
                plot!(
                    Dates.value.(ccrang)./1000,
                    wght_tmp,
                    label = "CC Weights",
                )
                plot!(
                    Dates.value.(ccrang)./1000,
                    tmp_ccs .* wght_tmp,
                    label = "Weighted CCs",
                )
            end
            hp_tmp_tr = plot(
                [xy_tmp,sac_tmp_segs[i_tmp]],
                label = ["xy" "sac"],
                xlabel = "samples",
            ) 
            hp_tmp_all = plot(hp_tmp_cc,hp_tmp_tr,layout=grid(2,1),size=(1200,1200))
            savefig(hp_tmp_all,string(cPlotOut,xy_name,"_cc_diag/",cc_wind_ctr[i],".pdf"))
        end
    else
        CC[i] = NaN
        Tdiff[i] = Tdiff[i-1]
    end
end
print("\n")

## LINEARLY INTERPOLATE TDIFF TO GET NEW TIMES
tmpt = vcat(
    0,
    Dates.value.(cc_wind_ctr .- xy_T[1]),
    Dates.value(xy_T[end]-xy_T[1]),
)
tmptdiff = vcat(0, Dates.value.(Tdiff), 0)
itp = LinearInterpolation(tmpt, tmptdiff)
xy_Tdiff = itp(Dates.value.(xy_T .- xy_T[1]))
xy_Tdiff = Dates.Millisecond.(round.(xy_Tdiff))
xy_Tcorr = xy_T .+ xy_Tdiff # correct xy  times
xy_Tdiff_prime = diff(Dates.value.(xy_Tdiff))./Dates.value.(diff(xy_T))

## GET TRACE ZERO LINES
print("Acquiring trace zero lines...\n")
trace_z = []
for i in ProgressBar(1:lastindex(trace_x))
    # get pseudo-linear trend (from start and ends of traces)
    local strt_end_idx = unique(vcat((1:dtrnd_wndw),(length(trace_x[i])-dtrnd_wndw:length(trace_x[i]))))
    strt_end_idx = strt_end_idx[findall(1 .<= strt_end_idx .<= length(trace_x[i]))]
    local gidx = intersect(
        findall(!isnan,trace_x[i][strt_end_idx]),
        findall(!isnan,trace_y[i][strt_end_idx])
    )
    local x = trace_x[i][strt_end_idx[gidx]]
    local y = trace_y[i][strt_end_idx[gidx]]
    local a = ((length(x)*sum(x.*y)-(sum(x)*sum(y))) /  # solve for slope
        ((length(x)*sum(x.^2))-(sum(x)^2)))
    local b = (1/(length(x)))*(sum(y)-a*sum(x)) # solve for y-intercept
    tmp_z = a*trace_x[i] .+ b # get trend
    if length(tmp_z)>(dtrnd_wndw*2)
        # get cosine-tapered moving-average
        tmp_y = trace_y[i] .- tmp_z # subtract trend to get moving avg.
        local dtrnd_strt = convert.(Int,1:dtrnd_step:(length(trace_x[i])-dtrnd_wndw)) # get starts (for computing values)
        local dtrnd_ctr = convert.(Int,round.(dtrnd_strt .+ dtrnd_wndw/2)) # get centers (for interpolating)
        local dtrnd_y = map(x->
            mean(filter(!isnan,tmp_y[dtrnd_strt[x]:(dtrnd_strt[x]+dtrnd_wndw)])),
            1:lastindex(dtrnd_strt),
        ) 
        local itp_y = LinearInterpolation(dtrnd_ctr, dtrnd_y) # get interpolant
        local trend_y = itp_y(dtrnd_ctr[1]:dtrnd_ctr[end]) # interpolate
        local trend_y_tapered = vcat( 
            zeros(dtrnd_ctr[1]-1),
            trend_y.*vcat( # cosine taper with flat top
                (cos.(range(pi,2*pi,length=dtrnd_ctr[1])).+1)./2, # from -1 to 1
                ones(length(trend_y)-2*dtrnd_ctr[1]),
                (cos.(range(0,pi,length=dtrnd_ctr[1])).+1)./2
            ), # from 1 to -1 
            zeros(length(tmp_y)-dtrnd_ctr[end]),
        ) 
        # add linear and moving average
        tmp_z = tmp_z .+ trend_y_tapered
    end
    # save the trace zero
    push!(trace_z, tmp_z)
end
print("\n")

## PLOT RESULTS
print("Plotting results...\n")
# plot trace time diff as seismogram !!! NOTE DOWN IS UP !!!
hp_tdiff = plot(
    trace_x,
    trace_y,
    la = 0.5,
    legend = false,
    yflip = true,
    title = string(xy_name," Tdiff"),
    size = (im_lgth,im_hght),
)
plot!(
    hp_tdiff,
    trace_x,
    trace_z,
    lc = :black,
    ls = :dot,
)
# sort Tdiff into lines
plt_Tdiff = []
plt_TdiffX = []
global idx_strt = 1
for i = 1:lastindex(trace_x)
    append!(plt_Tdiff, 
        trace_z[i].-Dates.value.(
            xy_Tdiff[idx_strt:(idx_strt+length(trace_x[i])-1)])./100)
    append!(plt_TdiffX, trace_x[i])
    append!(plt_Tdiff, NaN)
    append!(plt_TdiffX, NaN)
    global idx_strt = length(trace_x[i])+idx_strt
end 
# plot them
plot!(
    hp_tdiff,
    plt_TdiffX,
    plt_Tdiff,
    lc=:black,
)
savefig(hp_tdiff,string(cPlotOut,xy_name,"_sgram_tdiff.pdf"))

# crop sac
sac_Dcrp = sac_D[findall(xy_Tcorr[1] .<= sac_T .<= xy_Tcorr[end])]
sac_Tcrp = sac_T[findall(xy_Tcorr[1] .<= sac_T .<= xy_Tcorr[end])]
# detrend sac
sac_Dcrp = sac_Dcrp .- mean(filter(!isnan,sac_Dcrp))
# interpolate to length of xy
sac_range = range(0,Dates.value(sac_Tcrp[end]-sac_Tcrp[1]),length=length(xy_Y))
sac_Tcrp = sac_Tcrp[1] .+ Dates.Millisecond.(round.(sac_range))
sac_itp = LinearInterpolation(
    range(1,length(xy_Y),length=length(sac_Dcrp)),
    sac_Dcrp,
)
sac_Dcrp = sac_itp(1:length(xy_Y))
# plot time corrected xy along with sac
hp_trc = plot(
    xy_Tcorr,
    xy_Y./median(filter(!isnan,abs.(xy_Y))),
    label = "xy_corr",
    lc = :blue,
    title = xy_name,
    size = (1800,400),
)
plot!(
    hp_trc,
    xy_T,
    xy_Y./median(filter(!isnan,abs.(xy_Y))),
    lc = :blue,
    ls = :dot,
    label = "xy_raw",
)
plot!(
    hp_trc,
    sac_Tcrp,
    sac_Dcrp./median(filter(!isnan,abs.(sac_Dcrp))),
    lc = :red,
    label = "sac",
)
# plot time diff 1st deriv and cc on time axis
hp_diff = plot(
    xy_Tcorr,
    Dates.value.(xy_Tdiff)./1000,
    legend = false,
    title = "Tdiff in Seconds",
    size = (1800,400),
)
hp_fdt = plot(
    xy_Tcorr[1:end-1],
    xy_Tdiff_prime,
    legend = false,
    title = "First derivative (seconds/second)",
    size = (1800,400),
)
hp_cc = plot(
    cc_wind_ctr,
    CC,
    legend = false,
    title = "Cross Corr Coeff",
    size = (1800,400),
)
# aggregate
hp_all = plot(
    hp_trc,
    hp_diff,
    hp_fdt,
    hp_cc,
    layout = grid(4,1),
    size = (1800,1800)
)
savefig(hp_all,string(cPlotOut,xy_name,"_xy_sac.pdf"))
# make small seg plots
if plot_small_segs
    print("Plotting small segments...\n")
    if !isdir(string(cPlotOut,xy_name,"_xy_sac/"))
        mkdir(string(cPlotOut,xy_name,"_xy_sac/"))
    end
    small_seg_strt = xy_Tcorr[1]:small_seg_wndw:(xy_Tcorr[end]-small_seg_wndw)
    for i in ProgressBar(1:lastindex(small_seg_strt))
        # set bounds
        plot!(
            hp_trc, 
            xlims = (small_seg_strt[i],small_seg_strt[i]+small_seg_wndw),
            ylims = (
                minimum(
                    [minimum(filter(
                        !isnan,xy_Y[findall(
                            small_seg_strt[i] .< xy_Tcorr .< small_seg_strt[i]+small_seg_wndw)]
                    ))*1.5/median(filter(!isnan,abs.(xy_Y))),
                    minimum(filter(
                        !isnan,sac_D[findall(
                            small_seg_strt[i] .< sac_T .< small_seg_strt[i]+small_seg_wndw)]
                    ))*1.5/median(filter(!isnan,abs.(sac_D)))]
                ),
                maximum(
                    [maximum(filter(
                        !isnan,xy_Y[findall(
                            small_seg_strt[i] .< xy_Tcorr .< small_seg_strt[i]+small_seg_wndw)]
                    ))*1.5/median(filter(!isnan,abs.(xy_Y))),
                    maximum(filter(
                        !isnan,sac_D[findall(
                            small_seg_strt[i] .< sac_T .< small_seg_strt[i]+small_seg_wndw)]
                    ))*1.5/median(filter(!isnan,abs.(sac_D)))]
                )
            ),
        )
        plot!(
            hp_diff, 
            xlims = (small_seg_strt[i],small_seg_strt[i]+small_seg_wndw),
            ylims = (
                minimum(filter(
                    !isnan,(Dates.value.(xy_Tdiff)./1000)[findall(
                        small_seg_strt[i] .< xy_Tcorr .< small_seg_strt[i]+small_seg_wndw)]
                ))-0.5,
                maximum(filter(
                    !isnan,(Dates.value.(xy_Tdiff)./1000)[findall(
                        small_seg_strt[i] .< xy_Tcorr .< small_seg_strt[i]+small_seg_wndw)]
                ))+0.5
            ),
        )
        plot!(
            hp_fdt, 
            xlims = (small_seg_strt[i],small_seg_strt[i]+small_seg_wndw),
            ylims = (
                minimum(filter(
                    !isnan,xy_Tdiff_prime[findall(
                        small_seg_strt[i] .< xy_Tcorr[1:end-1] .< small_seg_strt[i]+small_seg_wndw)]
                ))-0.05,
                maximum(filter(
                    !isnan,xy_Tdiff_prime[findall(
                        small_seg_strt[i] .< xy_Tcorr[1:end-1] .< small_seg_strt[i]+small_seg_wndw)]
                ))+0.05
            ),
        )
        plot!(
            hp_cc, 
            xlims = (small_seg_strt[i],small_seg_strt[i]+small_seg_wndw),
            ylims = (
                minimum(filter(
                    !isnan,CC[findall(
                        small_seg_strt[i] .< cc_wind_ctr .< small_seg_strt[i]+small_seg_wndw)]
                ))-0.15,
                maximum(filter(
                    !isnan,CC[findall(
                        small_seg_strt[i] .< cc_wind_ctr .< small_seg_strt[i]+small_seg_wndw)]
                ))+0.15
            ),
        )
        # aggregate
        local hp_all = plot(
            hp_trc,
            hp_diff,
            hp_fdt,
            hp_cc,
            layout = grid(4,1),
            size = (1800,1800)
        )
        savefig(hp_all,string(cPlotOut,xy_name,"_xy_sac/",Dates.format(small_seg_strt[i],"yyyymmdd_HHMM"),".pdf"))
    end
end

# plot histogram of timediffs
hp_hdiff = histogram(
    Dates.value.(xy_Tdiff)./1000,
    nbins=150,
    title="Tdiff (sec)",
)
hp_hfdt = histogram(
    diff(Dates.value.(xy_Tdiff))./(Dates.value.(diff(xy_Tcorr))),
    nbins=150,
    title=string("First Deriv. (sec/sec) (mean = ",mean(abs.(xy_Tdiff_prime)),")"),
)
hp_hcc = histogram(CC,nbins=150,title="Cross Corr Coeff")
hp_hall = plot(hp_hdiff,hp_hfdt,hp_hcc,layout=grid(3,1),size=(1000,1200))
savefig(hp_hall,string(cPlotOut,xy_name,"_hists.pdf"))
print("\n")

## SAVE THE NEW TIMES
# save the times as milliseconds since start time
stime = minimum([xy_Tcorr[1],xy_T[1]])
xy_T_ms = Dates.value.(xy_T .- stime)
xy_T_ms_corr = Dates.value.(xy_Tcorr .- stime)
# print
f = open(string(cDataOut,xy_name,"_timed.txt"),"w")
print(f,string("stime = ",stime,"\n"))
print(f,"Times are given in milleconds past stime\n")
print(f,"Both linearly interpolated and cross-correlation corrected times are given\n")
print(f,string("stitched_t_linear,stitched_t_corr,stitched_x,stitched_y,detrend_y,orig_x,orig_y,orig_z,\n"))
global idx_strt = 1
for i = 1:lastindex(trace_x)
    for j = 1:lastindex(trace_x[i])
        tmpidx = idx_strt + j - 1
        print(f,string(
            xy_T_ms[tmpidx],",",
            xy_T_ms_corr[tmpidx],",",
            xy_X[tmpidx],",",
            xy_Y0[tmpidx],",",
            xy_Y[tmpidx],",",
            trace_x[i][j],",",
            trace_y[i][j],",",
            trace_z[i][j],",\n"
        ))
    end
    if i < length(trace_x)
        # update idx_strt
        global idx_strt = length(trace_x[i])+idx_strt
        # add NaN row
        print(f,string(
                xy_T_ms[idx_strt],",",
                xy_T_ms_corr[idx_strt],",",
                NaN,",",
                NaN,",",
                NaN,",",
                NaN,",",
                NaN,",",
                NaN,",\n"
            ))
    end
end 
close(f)
print(string("Saved detrended and stitiched xy_data to ",cDataOut,xy_name,"_timed.txt\n"))
print("\n")

## SAVE SCRIPT WITH PARAMETERS
f = open("/Users/thomaslee/Research/Microform/compare_xy_sac.jl","r")
ln = readlines(f)
close(f)
fset = open(string(cDataOut,xy_name,"_setup.txt"),"w")
global print_step = 0 # 0 is not printing, 1 is printing, 2 is done
global ln_count = 1
while print_step < 2
    if print_step == 0
        if occursin("## SETUP",ln[ln_count]) # check if we should start printing
            print(fset,string(ln[ln_count],"\n"))
            global print_step += 1 # make print_step 1
        end
    elseif print_step == 1
        if occursin("## SETUP DIRECTORIES",ln[ln_count]) # check if done printing
            global print_step += 1 # make print_step 2
        else # add line
            print(fset,string(ln[ln_count],"\n"))
        end
    end
    global ln_count += 1
end
close(fset)

print("Done!\n")
