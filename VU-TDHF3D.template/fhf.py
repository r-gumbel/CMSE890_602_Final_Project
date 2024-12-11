#!/usr/bin/env python3
import subprocess
f1m = 4
f1p = 2
f2m = 208
f2p = 82
directory = "12C24O"
for ri in range(200,30,-5):
    r=float(ri)/10
    f=open('tdhf3d.inp','w')
    inputfile ="""SLy4dL
 2                                             nof
 1                                             imode
  0.0    0.0                                   centf(1,if) if=1,nof
  0.0    0.0                                   centf(2,if) if=1,nof
  0.0    0.0                                   centf(3,if) if=1,nof
  0.0    0.0                                   boostf(1,if) if=1,nof
  0.0    0.0                                   boostf(2,if) if=1,nof
  0.0    0.0                                   boostf(3,if) if=1,nof
 00.0   00.0                                   euler_alpha(if)
 00.0   00.0                                   euler_beta(if)
 00.0   00.0                                   euler_gamma(if)
 """+str(f1m)+""".0D0 """+str(f2m)+""".0D0                                fmass(if) if=1,nof
 """+str(f1p)+""".0D0 """+str(f2p)+""".0D0                                 fcharg(if) if=1,nof
 3.2 3.4                                       radinf(1,if) if=1,nof
 3.2 3.1                                       radinf(2,if) if=1,nof
 3.2 3.2                                       radinf(3,if) if=1,nof
 0.020   35.0                                  x0dmp,e0dmp
  0000  2000  6.0D-3 2.0D-2                   itrbx,mtrbx,serr,derr
  1   20   00   1                               iprint,mprint,mplots,mplott
  0000 0 2000 100                              irest,irwgs,mrests,mrestt
  0 1                                          ifixcm,icoul
  0 0.000 1  1.90 5e-5                               iconstr,q2in,mconstr,c0=(1.9),d0=(5e-5)
  0 0 1 0                                      itimrev, itimrevs, iodds, hdiag
  0 0                                          ipairf(1,if) if=1,nof
  0 0                                          ipairf(2,if) if=1,nof
  0 2000                                       nexadd, nexiter
  26 20                                         nextra_n(if), if=1,nof
  00 20                                         nextra_p(if), if=1,nof
 1                                             ifixb
 00.0  """+str(r)+"""   0.0                            ecm,rsep,xb
 12     1.0D-17                                mxp,terr
 1      0.400D0                                nt,dt"""
    f.write(inputfile)
    f.close()
    subprocess.call("./run", shell=True)
    subprocess.call("cp -r AL_"+str(f1m)+"_ZL_"+str(f1p)+"_AR_"+str(f2m)+"_ZR_"+str(f2p)+"_SLy4dL_Ecm0.000_beta0_b0.000" \
            " AL_"+str(f1m)+"_ZL_"+str(f1p)+"_AR_"+str(f2m)+"_ZR_"+str(f2p)+"_SLy4dL_Ecm0.000_beta0_b0.000_R_"+str(r), shell=True)
