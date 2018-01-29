;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; R. Nair IMAGe/NCAR        log-plot for convergence                    ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;
;  Main program
;

begin

   m = 3  ;; For 3-pt plot



   xx     = new((/m/),float)
   f1     = new((/m/),float)
   f2     = new((/m/),float)
   o3     = new((/m/),float)
   xp     = new((/m/),float)

   n = 6  ;; For 3-pt plot
   x1     = new((/n/),float)
   xm     = new((/n/),float)
   y1     = new((/n/),float)
   y2     = new((/n/),float)
   o4     = new((/n/),float)

;; 3 pt convergence test
   xx(0) = 1
   xx(1) = 2
   xx(2) = 3
 ; xx(3) = 4

   xp(0) = 4
   xp(1) = 6
   xp(2) = 8
;  xp(3) = 10

;; SD L1 data
    f1(0) = 1.290760E-03
    f1(1) = 4.290999E-05
    f1(2) = 3.854418E-06
 ;  f1(3) = 3.667126E-12

;; SD linf data
    f2(0) = 1.296129E-03
    f2(1) = 2.946655E-05
    f2(2) = 7.022191E-06
 ;  f2(3) = 9.213971E-11

     t=5.0E-02
;  o4(0) = t/(5.)^3
;  o4(1) = t/(7.)^3
;  o4(2) = t/(9.)^3
;  o4(3) = t/(11.)^3
;  o4(4) = t/(13.)^3
;  o4(5) = t/(15.)^3

   dat2     = new((/2,m/),float)
   dat2(0,:) = f1
   dat2(1,:) = f2

;  dat2(1,:) = o3
;  dat2(2,:) = o4

;; For second plot on the right

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; create plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;

  X11=1      ;; for X-window
  esp=0      ;; for eps file
  PDF=0

  plots = new(1,graphic)

  if (X11 .eq. 1) then
    wks = gsn_open_wks("x11","hp_sw2_sd")
  end if
  if (esp .eq. 1) then
    wks = gsn_open_wks("eps","hp_sw2_sd")
  end if
  if (PDF .eq. 1) then
    wks = gsn_open_wks("pdf","hp_sw2_sd")
  end if

  res               = True                   ; plot mods desired
  res@gsnDraw     = True      ; Will draw later, after overlaying
  res@gsnFrame    = False      ; all plots

  res@tiMainString  = "Gaussian Hill"   ; title
  res@tiYAxisString = "Max Value"           ; y axis title
  res@tiXAxisString = "Time"                 ; x axis title

  res@trXMinF           = 0.5
  res@trXMaxF           = 3.5
  res@trYMinF           = 1.0E-07
  res@trYMaxF           = 1.0E-02

  res@tmXBMode   = "Explicit"
  res@tmXBValues = xx                       ; location of labels ^M
  res@tmXBLabels = xp
  res@tmLabelAutoStride = True                    ; nice stride on labels
  res@xyLineLabelFontAspectF =2.0

 res@tiXAxisString = "Nv-1(Poly. Degree)"


 ;;Change the Y axis scale to be log and redraw the plot.
  res@trYLog        = True
; res@tiMainString  = "(a) L~B1~2  Error"
  res@tiMainString  = "NH-IGW 3D Test"
  res@tiYAxisString = "Normalized Errors: L~B1~~F18~x"
  res@xyComputeYMin = True   ; To keep from getting some warning msgs.


  res@xyLineColors      = (/"black","black"/)  ; line colors
  res@xyLineThicknesses = (/2.0,2.0/)        ; line thicknesses
  res@xyDashPatterns    = (/0.0,1.0,2.0/)        ; line patterns
  res@xyMarkLineModes   = (/"MarkLines","MarkLines"/)
; res@xyMarkLineModes   = (/"MarkLines"/)
  res@xyMarkerColors    = (/"black","black"/)   ;
  res@xyMarkers         = (/16,5/)          ; style
  res@xyMarkerSizeF = 0.0125


;; Legend markers

  res@pmLegendDisplayMode    = "Always"            ; turn on legend^M
  res@pmLegendSide           = "Top"               ; Change location of ^M
  res@pmLegendParallelPosF   = .30                  ; move units right^M
  res@pmLegendOrthogonalPosF = -1.1               ; move units down^M
  res@pmLegendWidthF         = 0.30                ; Change width and^M
  res@pmLegendHeightF        = 0.10                ; height of legend.^M
  res@lgPerimOn              = False               ; turn off box around^M
  res@lgLabelFontHeightF     = 16                ; label font height

  res@lgPerimFill            = "SolidFill"      ; Fill legend box w/white
  res@lgPerimFillColor       = "white"          ; so it masks XY curves

  res@xyLineLabelFontHeightF = 15.0        ; Font size and color
  res@xyLineLabelFontColor   = 10           ; for line labels
  res@xyLineLabelFontAspectF =2.0
; res@xyExplicitLegendLabels  = (/" W5"," KL","O(3)","O(4)"/)         ; create explicit labels
; res@xyExplicitLegendLabels  = (/" W5","O(2)","O(3)"/)         ; create explicit labels
  res@xyExplicitLegendLabels  = (/"L~B1~1","L~B1~~F18~x"/)         ; create explicit labels

 ;plots(0) = gsn_csm_xy(wks,xx,dat2,res)          ; Multiplot
  plot = gsn_csm_xy(wks,xx,dat2,res)          ; Multiplot


frame(wks)

end


