;-----------------------------------------------------------------------------------------
; NAME:     orthreg_lgcens                                          IDL Function
;
; PURPOSE: Model function for an orthogonal 
;
; CALLING SEQUENCE: 
; parms = MPFIT('orthreg_lgcens', p0, $
;             FUNCTARGS={X: X, Y: Y, SIGMA_X: SIGMA_X, SIGMA_Y: SIGMA_Y, rho: rho_err}, $
;             ...)
;
; INPUTS: 
;     p0 - array containing an initial guess for the number of fitted parameters
;           (here three: intercept, slope, intrinsic scatter)
;     Functargs - brings in the data values, errors, and flags as well as rho (a variable
;       containing the correlation of measurement errors)
;
; OUTPUTS: returns the inverse likelihoods (i.e., the value to be minimized by MPFIT)
;
; COMMENTS:
; 
; Note that the user is not meant to call this function directly, but instead pass
; it on to MPFIT
; 
; written based on the example of linfitex.pro (part of the mpfit package, available at
;         https://www.physics.wisc.edu/~craigm/idl/fitting.html)
; using the likelihood formula of Pihajoki 2017 arxiv:1704.05466v2, specifically
; equations 33-35 and B5-B6 from the Appendix.  
;
; REVISION HISTORY:
;   2018-Oct Written by Lauranne Lanz (Dartmouth)
;-----------------------------------------------------------------------------------------
function orthreg_lgcens, p, $
  x=x, y=y, sigma_x=sigma_x, sigma_y=sigma_y, flg_x=flg_x, flg_y=flg_y,rho = rho_err, $
  _EXTRA=extra
  
  a = p[0]   ;; Intercept
  b = p[1]   ;; Slope
  ins = p[2]  ;; intrinsic scatter

  ndt = n_elements(x)
  lklihd = make_array(ndt, /double)


  for i = 0, ndt-1 do begin

    ; figure out which likelihood function to use
    flg = 0
    if flg_x[i] eq 0 then begin
      if flg_y[i] eq 0 then flg = 0 else flg = 1
    endif else begin
      if flg_y[i] eq 0 then flg = 2 else flg = 3
    endelse

    case flg of
      0: begin ; both directions detected -> use eqn 33-35
        nu2_i = (y[i] - a - b*x[i])^2 / (1 + b^2) ; eqn 34
        rho_i = rho_err
        sig2_i = (b^2 * sigma_x[i]^2 + sigma_y[i]^2 - 2*b*rho_i*sigma_x[i]*sigma_y[i])/(1 + b^2) + ins^2 ; eqn 35
        lklihd[i] = exp(-nu2_i/sig2_i) /sqrt(2*!pi*sig2_i) ; eqn 33
      end
      1: begin ; x is detected but y is a UL use equation B6
        bk = 10. ; the logarithmic base of our measurements
        term1 = a + x[i]*b + 1/2*(b^2*sigma_x[i]^2 + (1+b^2)*ins^2)*alog10(bk)
        term2 = erfc( (a + x[i]*b + (sigma_y[i]^2 + (1+b^2)*ins^2)*alog10(bk) - y[i] )/(sqrt(2 * (b^2*sigma_x[i]^2 + (1+b^2)*ins^2))) )
        lklihd[i] = alog10(bk)*sqrt((1+b^2))/(2*y[i])* bk^(term1) * term2
        if lklihd[i] eq 0 then lklihd[i]=1e-5       ; prevent an issue when numbers get too small  
      end
      2: begin ; y is detected but x is a UL use equation B5
        bk = 10. ; the logarithmic base of our measurements
        term1 = (2 * (y[i] - a)*b + (sigma_y[i]^2 + (1+b^2)*ins^2)*alog10(bk))/(2*b^2)
        term2 = erfc( ((y[i]-a)*b + (sigma_y[i]^2 + (1+b^2)*ins^2)*alog10(bk) - b^2* x[i] )/(sqrt(2*(sigma_y[i]^2 + (1+b^2)*ins^2))*abs(b)) )
        lklihd[i] = alog10(bk)*sqrt((1+b^2)/b^2)/(2*x[i])* bk^(term1) * term2
        if lklihd[i] eq 0 then lklihd[i]=1e-5    ; prevent an issue when numbers get too small
        ;stop
      end
      ; Pihajoki did not include a formula for when both X and Y are limits, so we effectively exclude these data points
      3: begin ; neither y nor x is detected, make the likelihood very high so that inverse essentially excludes it
        lklihd[i] = 1d10
      end
    endcase

  endfor

  invlklihd = 1./lklihd

  return, invlklihd
end
