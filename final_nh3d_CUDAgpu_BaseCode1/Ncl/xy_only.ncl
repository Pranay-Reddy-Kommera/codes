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
  fdat = "./data/data_Zcut.dat" ; plot XY-plane data

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

; Pick filetype from: pdf, ps, png, x11, ...
; filetype = "x11"
; wks = gsn_open_wks(filetype, "plot_rhoq")
  wks = gsn_open_wks("x11","plot_rhoq")

; gsn_define_colormap(wks,"rainbow+white")
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")

; Set up resources.
  res                     = True
  res@sfYArray              = zz
  res@sfXArray              = xx

  res@gsnMaximize         = True        ; Maximize plot in frame
  res@gsnSpreadColors     = True        ; Span full color map
  res@cnFillOn            = True        ; Turn on contour fill
  res@cnLinesOn           = True       ; Turn off/on contour lines
  res@lbLabelAutoStride   = True        ; Control labelbar spacing

  res@gsnSpreadColorStart = 10
  res@gsnSpreadColorEnd   = -4

  res@tiMainString          = "DG 3D advection"

  res@cnLevelSelectionMode = "ManualLevels"  ; manual contour levels
  res@cnMinLevelValF       = 0.1             ; contour levels.
  res@cnMaxLevelValF       = 1.0
  res@cnLevelSpacingF      = 0.1

 plot = gsn_csm_contour(wks,data,res)

end
