;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; R. Nair IMAGe/NCAR  Contour Advection Data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; ASCII data read
  gdat = "./data/grid_geometry.dat"
  xdat = "./data/grid_xs.dat"
  ydat = "./data/grid_ys.dat"
  zdat = "./data/grid_zs.dat"
  fdat = "./data/data_Ycut.dat" ; plot XZ-plane data
  mdat = "./data/mtn_surface.dat" ; plot XY-plane data

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
  m_dat = asciiread(mdat,(/ly,lx/),"float")

  xslice = m_dat(ly/2,:) / 1000. 

; Pick filetype from: pdf, ps, png, x11, ...
; filetype = "x11"
; wks = gsn_open_wks(filetype, "plot_rhoq")
  wks = gsn_open_wks("x11","plot_rhoq")

; gsn_define_colormap(wks,"rainbow+white")
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")

; Set up resources.
  res                     = True
  res@sfYArray              = zz /1000.
  res@sfXArray              = xx /1000.

  res@gsnMaximize         = True        ; Maximize plot in frame
  res@gsnSpreadColors     = True        ; Span full color map
  res@cnFillOn            = True        ; Turn on contour fill
  res@cnLinesOn           = True       ; Turn off/on contour lines
  res@lbLabelAutoStride   = True        ; Control labelbar spacing


  res@gsnDraw          = False
  res@gsnFrame         = False

  res@gsnSpreadColorStart = 0
  res@gsnSpreadColorEnd   = -4
;  res@tmXBLabelFont = 25         ; Change the font.
;  res@tmXBLabelFontHeightF = 0.020 ; Change the font.

  res@vpXF      = 0.15   ; Change the size and location of the
  res@vpYF      = 0.85   ; plot on the viewport.
  res@vpWidthF  = .80          ; choose width of plot
  res@vpHeightF = .40          ; choose height of plot

; res@tiMainString          = "DG 3D advection"

  res@cnLevelSelectionMode = "ManualLevels"  ; manual contour levels
  res@cnMinLevelValF       = 0.1             ; contour levels.
  res@cnMaxLevelValF       = 1.0
  res@cnLevelSpacingF      = 0.1

; res@gsnTickMarksOn = True
; res@tmXBLabelStride = 2
; res@tiXAxisString = "Km"                   ; x-axis label
  res@tiYAxisString = "Km"              ; y-axis label

   res@tmYLOn        = True       ; Turn on right axis tick marks.
   res@tmYLMinorOn   = True  
   res@tmYRMinorOn   = True 
   res@tmXTMinorOn   = True      ; No minor tick marks.
   res@tmXBMinorOn   = True     ; No minor tick marks.

  res@lbLabelStride       = 1       ;
  res@lbTitleFontThicknessF   = .5
  res@lbLabelFontHeightF    = .016

 plot = gsn_csm_contour(wks,data,res)

  lres           = True  
  lres@gsnDraw          = False
  lres@gsnFrame         = False
  lres@vpXF      = 0.15   ; Change the size and location of the
  lres@vpYF      = 0.85   ; plot on the viewport.
  lres@vpWidthF  = .80          ; choose width of plot
  lres@vpHeightF = .40          ; choose height of plot
  lres@trYMinF = 0.0                ; Set min/max of X axis.
  lres@trYMaxF = 30.
  lres@trXMinF = 0.0                ; Set min/max of X axis.
  lres@trXMaxF = 120.
  lres@xyLineThicknesses	 = (/2./) 
  lres@xyLineColor	 = "blue" 

;  lres@tmYLMinorOn   = False      ; No minor tick marks.
;  lres@tmXTMinorOn   = False      ; No minor tick marks.
;  lres@tmXBMinorOn   = False     ; No minor tick marks.

 lres@gsnYRefLineColor      = (/"black"/)
 lres@gsnAboveYRefLineColor = 100   ;(/"black"/)
 lres@gsnYRefLine           = (/0.0 /)


  xline = xx/1000. 
  pxy  = gsn_csm_xy(wks,xline,xslice,lres) 

 overlay(plot,pxy)
 draw(plot) 


 txres = True
 txres@txFontHeightF = 0.024
 txres@txFontColor = 1
 txres@txJust = "CenterCenter"
 gsn_text_ndc(wks, "DG3D Advection: Flow over a Mountain (X-Z)", 0.55, 0.85, txres)

 frame(wks)


end
