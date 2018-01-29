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
  fdat = "../data/mtn_surface.dat" ; plot XY-plane data

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
  data = asciiread(fdat,(/ly,lx/),"float")

; Pick filetype from: pdf, ps, png, x11, ...
; filetype = "x11"
; wks = gsn_open_wks(filetype, "plot_rhoq")
  wks = gsn_open_wks("x11","plot_rhoq")

; gsn_define_colormap(wks,"rainbow+white")
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")

; Set up resources.
  res                     = True
  res@sfYArray              = yy/1000.
  res@sfXArray              = xx/1000.

  res@gsnMaximize         = True        ; Maximize plot in frame
  res@gsnSpreadColors     = True        ; Span full color map
  res@cnFillOn            = True        ; Turn on contour fill
  res@cnLinesOn           = True       ; Turn off/on contour lines
  res@lbLabelAutoStride   = True        ; Control labelbar spacing
  res@gsnFrame            = False

  res@vpXF      = 0.15   ; Change the size and location of the
  res@vpYF      = 0.85   ; plot on the viewport.
  res@vpWidthF  = .80          ; choose width of plot
  res@vpHeightF = .40          ; choose height of plot

  res@gsnSpreadColorStart = 2 
  res@gsnSpreadColorEnd   = -2

; res@tiMainString          = "Mountain Surface (witch of Agnesi)"

  res@cnLevelSelectionMode = "ManualLevels"  ; manual contour levels
  res@cnMinLevelValF       = 0.0             ; contour levels.
  res@cnMaxLevelValF       = 5000.0
  res@cnLevelSpacingF      = 500.0 
; res@cnMinLevelValF       = 0.0             ; contour levels.
; res@cnMaxLevelValF       = 0.6
; res@cnLevelSpacingF      = 0.05 

 plot = gsn_csm_contour(wks,data,res)

 txres = True
 txres@txFontHeightF = 0.024
 txres@txFontColor = 1
 txres@txJust = "CenterCenter"
 gsn_text_ndc(wks, "Mountain Surface (wich of Agnesi)", 0.5, 0.85, txres)

 frame(wks)

end
