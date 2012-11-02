pro plot_cov, ips, fileIn, nbins, scaleMin, scaleMax, title, coeff_flag

;Reads fileIn
cov = dblarr(nbins,nbins)
coeff = dblarr(nbins,nbins)
openr, lun, fileIn, /get_lun
readf, lun, cov


;Settings
if (ips) then begin
    psfile = "graph.ps"
    set_plot, 'ps'
    device, filename=psfile, bits_per_pixel=8,/color $
       ,xsize=14,ysize=16,yoffset=7

endif else begin
    set_plot, 'x'
    window, xsize=500,ysize=500, retain=2
endelse

loadct, 0
ncolors=256
colortable=33


;Corr coeff
for i=0,nbins-1 do begin
    for j=0,nbins-1 do begin	
        coeff(i,j) = cov(i,j)/sqrt(cov(i,i)*cov(j,j))
    endfor
endfor

if(coeff_flag) then begin
    cov = coeff
endif

;max = MAX(cov)
;min = MIN(cov)

max = 1.0
min = -1.0


;Rescaling
cov = (ncolors-1)*(cov-min)/(max-min)

;Color bar
bar =  bindgen(ncolors-1)#replicate(1B,10)

;Top = zp zs comparison and bottom = legend
!P.MULTI = [0,1,2]

;Upper panel---------------------------------------------------------

;General settings
!P.charsize  = 1.2
!y.charsize  = !P.charsize
!x.charsize  = !P.charsize
symsize      = 1.5

!P.thick     = 5
!x.thick     = !P.thick
!y.thick     = !P.thick
!P.charthick = !P.thick

plot,[0],[0],/nodata,xstyle=1+4,ystyle=1+4,POSITION=[0.15, 0.30, 0.90, 0.95],xrange=[scaleMin,scaleMax],yrange=[scaleMin,scaleMax],title=title
tv, cov, !x.window(0),!y.window(0),xsize=!x.window(1)-!x.window(0),ysize=!y.window(1)-!y.window(0),/normal

axis, XAXIS=0,xtitle=textoidl('\theta(deg)'),xstyle=1,xrange=[scaleMin,scaleMax],/xlog
axis, XAXIS=1,xtickformat="(A1)",xstyle=1,xrange=[scaleMin,scaleMax],/xlog
axis, YAXIS=0,ytitle="",ystyle=1,yrange=[scaleMin,scaleMax],/ylog
axis, YAXIS=1,ytickformat="(A1)",ystyle=1,yrange=[scaleMin,scaleMax],/ylog

;legend
plot,[0],[0],/nodata,/noerase $
,xstyle=1+4,ystyle=1+4,xrange=[min,max],yrange=[0.0,1.0] $
,position=[0.15,0.12,0.90,0.17]

tv, bar, !x.window(0),!y.window(0),xsize=!x.window(1)-!x.window(0),ysize=!y.window(1)-!y.window(0), /normal

axis, XAXIS=0,xstyle=1;,xticks=4;,xticklen=0.2
plots, [min,min,max,max],[0.0,1.0,1.0,0.0]

if (ips) then begin
    device, /close
    set_plot, 'x'
endif

end
