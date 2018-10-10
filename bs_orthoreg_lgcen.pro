;-----------------------------------------------------------------------------------------
; NAME:    bs_orthoReg_lgcen                                    IDL Procedure
;
; PURPOSE: Run a bootstrap analysis to determine the confidence range for
;         an orthogonally regressed linear fit to logarithmic values (e.g., luminosities) which have
;         have censored values in one or both parameters
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;    dfile   - name of ASCII file containing columns of X parameter, X parameter 
;             uncertainty, X parameter flag (0 = value, 1 = upper limit),
;             followed by the same for Y
;    nboots  - integer number of the number of bootstrap samples 
;    nme_out - name for this dataset 
;
; OPTIONAL INPUTS:
;    pth    - path to directory in which to create working directories and read in 
;                       dfile (assumed to be . if not specified)
;    flgow  - flag about overwriting output if already exists (0 = no, 1 = yes), 
;                       default to 0
;    rhoE   - correlation between measurement errors (see comments)
;    flgpl  - flag about whether to plot the results (0 = no, 1 = yes)
;                       default to 0
;
; OUTPUTS: IDL savefile containing the results of each fit trial as well as fit_stats,
;          the original data, and the random indices drawn to create the bootstrapped
;          data
;          
;          fit_stats (also printed to the screen) contains:
;           - orig. slope, median of bs slopes, stdev of bs slopes, 
;                 median mpfit unc. of slopes
;           - orig intercept, median of bs intercepts, stdev of bs intercepts, 
;                 median mpfit unc. of intercepts
;           - orig. scatter, median of bs scatter, stdev of bs scatter, 
;                 median mpfit unc. of scatter
; COMMENTS:
;
;  Orthogonal regression fitting differs from least squares fitting in that (1)
;  errors in both X and Y are taken into account and the fitted line is selected
;  to minimize the perpendicular distance between the data and line rather than
;  the distance in the Y direction only.
;
; written based on linfitex.pro (part of the mpfit package, available at 
;         https://www.physics.wisc.edu/~craigm/idl/fitting.html)
; using the likelihood formula of Pihajoki 2017 arxiv:1704.05466v2, specifically 
; equations 33-35 and B5-B6 from the Appendix.  We assume rho (correlation between 
; measurement errors) eq 0, since the luminosities used were measured independently.
;
;
; EXAMPLES: bs_orthoReg_lgcen, 'test_data.dat', 1000, 'test'
;
; PROCEDURES CALLED:
;   - mpfit programs
;   - ortho_lgcens: function that returns the likelihood using the Pihajoki formulae
;
; REVISION HISTORY:
;   2018-Oct  Written by Lauranne Lanz (Dartmouth)
;-----------------------------------------------------------------------------------------
pro bs_orthoReg_lgcen, dfile, nboots, nme_out, $
                        pth = pthin, flgow = flg_ow, flgpl = flg_pl, rhoE = rho_err, $
                        VERBOSE = verbose

  if keyword_set(verbose) then begin
     print, ' '
     print, 'begin bs_orthoReg_lgcen'
     tic
  endif

  if n_elements(pthin) eq 0 then pthin = './' ; default to current directory
  if n_elements(flg_ow) eq 0 then flg_ow = 0  ; default to stop if going to overwrite
  if n_elements(flg_pl) eq 0 then flg_pl = 0  ; default to not plotting the results
  if n_elements(rho_err) eq 0 then rho_err = 0  ; for our purposes, we assumed the measurement errors were uncorrelated
  nboots = fix(nboots) ; make nboots into integer if it wasn't already
  if nboots gt 10000 then begin
    print, 'This many bootstrap samples may cause some unintended crashes. Check the results carefully.'
    print, 'Try no more than 1000 to get an initial sense.'
    stop
  endif
  nme_out = '_'+nme_out

  ndt = 7 ; slope, unc slope, intercept, unc intercept, intrinsic scatter, unc scatter, chi2/dof
  fit_bs = make_array(nboots+1, ndt)

  ; Check whether we need to create a working directory
  out_dir = pthin + 'bs_oreg'+nme_out+'/'
  if file_exists(out_dir) eq 0 then begin ; doesn't exist yet
    spawn, 'mkdir ' + out_dir
  endif else begin
    if flg_ow eq 0 then begin
      print, 'Warning: An output directory with this name already exists.'
      stop
    endif
  endelse

  ; import data 
  readcol, pthin+dfile, val_X, err_X, fg_X, val_Y, err_Y, fg_Y, format='F,F,I,F,F,I', /silent
  n_mt = n_elements(val_X)

  iboots = round(randomu(sd, n_mt, nboots)*(n_mt-1))

  ; loop over the bootstrap samples
  for bs=0, nboots do begin

    if bs eq 0 then begin ; this is the original one
      ft_i = mpfit('orthreg_lgcens', [1d, 1d, 1d], functargs={x: val_X, y: val_Y, sigma_x: err_X, sigma_y: err_Y, $
                    flg_x: fg_X, flg_y: fg_Y, rho: rho_err}, $
                    perror = ft_iu, bestnorm = ft_i_chi, dof=ft_i_dof, /quiet)
    endif else begin
      bs2 = bs - 1 ; to avoid forgetting the -1 somewhere
      rand_ind = reform(iboots[*,bs2])
      X_bs  = val_X[rand_ind]
      err_X_bs = err_X[rand_ind]
      fg_X_bs = fg_X[rand_ind]
      val_Y_bs  =val_Y[rand_ind]
      err_Y_bs = err_Y[rand_ind]
      fg_Y_bs = fg_Y[rand_ind]      
      ft_i = mpfit('orthreg_lgcens', [1d, 1d, 1d], functargs={x: X_bs, y: val_Y_bs, sigma_x: err_X_bs, sigma_y: err_Y_bs, $
                    flg_x: fg_X_bs, flg_y: fg_Y_bs, rho: rho_err}, $
                    perror = ft_iu, bestnorm = ft_i_chi, dof=ft_i_dof, /quiet)
    endelse

    ; store the mpfit outputs
    fit_bs[bs, 0] = ft_i[1]  ; Slope
    fit_bs[bs, 1] = ft_iu[1]  ; Unc Slope
    fit_bs[bs, 2] = ft_i[0]  ; Intercept
    fit_bs[bs, 3] = ft_iu[0]  ; Unc Intercept
    fit_bs[bs, 4] = ft_i[2]  ; Intrinsic Scatter
    fit_bs[bs, 5] = ft_iu[2]  ; Unc. Intrinsic Scatter
    fit_bs[bs, 6] = ft_i_chi/ft_i_dof  ; inv. likelihood ^2 / dof
      
    if keyword_set(verbose) and (bs mod 10 eq 0) then begin
      print, 'Reached bootstrap sample: ', bs
      
      toc
    endif
  endfor 
   
  ; orig. slope, median of bs slopes, stdev of bs slopes, median mpfit unc. of slopes
  ; orig intercept, median of bs intercepts, stdev of bs intercepts, median mpfit unc. of intercepts
  ; orig. scatter, median of bs scatter, stdev of bs scatter, median mpfit unc. of scatter
  fit_stats = [[fit_bs[0,0], median(fit_bs[1:*, 0]), stdev(fit_bs[*, 0]), median(fit_bs[1:*, 1])], $
               [fit_bs[0,2], median(fit_bs[1:*, 2]), stdev(fit_bs[1:*, 2]), median(fit_bs[1:*, 3])], $
               [fit_bs[0,4], median(fit_bs[1:*, 4]), stdev(fit_bs[1:*, 4]), median(fit_bs[1:*, 5])]]

  save, filename=out_dir+'ortho_reg_'+nme_out+'.sav', fit_bs, fit_stats, iboots, val_X, err_X, fg_X, val_Y, err_Y, fg_Y
  
  print, 'Results: '
  print, 'Slope    : '
  print, 'Original Data: ', fit_stats[0, 0], '   Median of BS trials: ', fit_stats[1, 0]
  print, 'Stdev of BS trials: ', fit_stats[2, 0], '   Median of Mpfit Slope unc: ', fit_stats[3, 0]
  print, ' '
  print, 'Intercept    : '
  print, 'Original Data: ', fit_stats[0, 1], '   Median of BS trials: ', fit_stats[1, 1]
  print, 'Stdev of BS trials: ', fit_stats[2, 1], '   Median of Mpfit Intercept unc: ', fit_stats[3, 1]
  print, ' '
  print, 'Intrinsic Scatter    : '
  print, 'Original Data: ', fit_stats[0, 2], '   Median of BS trials: ', fit_stats[1, 2]
  print, 'Stdev of BS trials: ', fit_stats[2, 2], '   Median of Mpfit Scatter unc: ', fit_stats[3, 2]

  if flg_pl eq 1 then begin
    set_plot, 'x'
    plot, val_X, val_Y, psym=4, xr=[min(val_X)-0.5, max(val_X)+0.5], yr=[min(val_Y)-0.5, max(val_Y)+0.5]
    srtx = val_X[sort(val_X)]
    for bs = 1, nboots do oplot, srtx, fit_bs[bs, 0]*srtx + fit_bs[bs, 2], linestyle = 0, thick=0.5
    oplot, srtx, median(fit_bs[1:*, 0])*srtx + median(fit_bs[1:*, 2]), linestyle = 0, thick=5, color=cgcolor('red')
    oplot, srtx, fit_bs[0, 0]*srtx + fit_bs[0, 2], linestyle = 0, thick=5, color=cgcolor('blue')
    oplot, val_X, val_Y, psym=4, symsize=1, color=cgcolor('cyan'), thick=2
  endif   
    
  if keyword_set(verbose) then begin
    print, 'Reached end of bs_orthoReg_lgcen '
    toc
  endif
 stop
end