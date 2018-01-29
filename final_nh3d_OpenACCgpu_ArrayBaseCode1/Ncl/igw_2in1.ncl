;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; R. Nair IMAGe/NCAR  Contour Advection Data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; ASCII data read
  gdat = "../data/grid_geometry.dat"
  xdat = "../data/grid_xs.dat"
  ydat = "../data/grid_ys.dat"
  zdat = "../data/grid_zs.dat"
  fdat = "../data/data_Ycut.dat" ; plot XZ-plane data

; Dimension of the data (lat/lon)
  gd = asciiread(gdat,(/4/),"integer")
  nx = gd(0)
  nex = gd(1)
  ney = gd(2)
  nez = gd(3)

  lx = (nx-1)*nex + 1
  ly = (nx-1)*ney + 1
  lz = (nx-1)*nez + 1

  xx = asciiread(xdat,(/lx/),"float")
  yy = asciiread(ydat,(/ly/),"float")
  zz = asciiread(zdat,(/lz/),"float")
  data = asciiread(fdat,(/lz,lx/),"float")
 ff     = new((/lx/),float)

 ff(:) = data(lz/2,:)    ;z-half

; Pick filetype from: pdf, ps, png, x11, ...
; filetype = "x11"
; wks = gsn_open_wks(filetype, "plot_rhoq")
  wks = gsn_open_wks("x11","plot_rhoq")

; gsn_define_colormap(wks,"rainbow+white")
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")

  res0 = True   ;; common 'resources'

  res0@gsnFrame           = False   ; Don't advance the frame
  res0@vpHeightF          = 0.20    ; Width and height of
  res0@vpWidthF           = 0.75    ; each plot in NDC space.

  res = res0

; Set up resources.
  res                     = True
  res@sfYArray              = zz /1000.
  res@sfXArray              = xx /1000.

  res@gsnMaximize         = False        ; Maximize plot in frame
  res@gsnSpreadColors     = True        ; Span full color map
  res@cnFillOn            = True        ; Turn on contour fill
  res@cnLinesOn           = True       ; Turn off/on contour lines
  res@lbLabelAutoStride   = True        ; Control labelbar spacing

  res@gsnFrame = False

  res@gsnSpreadColorStart = 35
  res@gsnSpreadColorEnd   = -1
;  res@tmXBLabelFont = 25         ; Change the font.
;  res@tmXBLabelFontHeightF = 0.020 ; Change the font.

  res@vpXF      = 0.15   ; Change the size and location of the
  res@vpYF      = 0.85   ; plot on the viewport.
; res@vpWidthF  = .75          ; choose width of plot
; res@vpHeightF = .30          ; choose height of plot

; res@tiMainString          = "DG 3D advection"

  res@cnLevelSelectionMode = "ManualLevels"  ; manual contour levels
  res@cnMinLevelValF       = -0.005             ; contour levels.
  res@cnMaxLevelValF       = 0.005
  res@cnLevelSpacingF      = 0.001 

  res@lbLabelStride       = 1       ;
  res@lbTitleFontThicknessF   = .5
  res@lbLabelFontHeightF    = .016
  res@tiYAxisString = "z (Km)"

 plots = gsn_csm_contour(wks,data,res)

;; Lower plot
;;=========================
  res1 = True
  res1 = res0

; res1@gsnFrame         = False
  res1@xyLineColor       = "blue"
  res1@trYMinF =-0.0060                ; Set min/max of X axis.
  res1@trYMaxF = 0.0060
  res1@tiMainString  = "Potential Temp Perturbation (~F33~q~F33~'~F21~)  at z=5km"
  res1@tiMainFontHeightF = 0.022
  res1@tiXAxisString = "x (km)"
  res1@xyLineThicknesses = 2.0

 res1@vpXF      = 0.15  ; Change the size and location of the
 res1@vpYF      = 0.4   ; plot on the viewport.

  xs = xx/ 1000.0 
 plots = gsn_csm_xy(wks,xs,ff,res1)

 txres = True
 txres@txFontHeightF = 0.018
 txres@txFontColor = 1
 txres@txJust = "CenterCenter"
 gsn_text_ndc(wks, " HOMAM-C (F3D): IGW (~F33~q~F33~'~F21~), [np=4,Elem=64x32x20,t=3000s]", 0.5, 0.88, txres)

 frame(wks)


end
