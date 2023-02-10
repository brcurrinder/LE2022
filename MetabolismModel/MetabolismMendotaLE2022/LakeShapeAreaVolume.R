LakeShapeAreaVolume = function(Zmax,Zmean,Ztherm,Plot){
  # Algorith for lake shape
  # from: Carpenter SR. 1983. Lake geometry: Implications for production and sediment accretion rates. 
  #            J. Theor. Biol. 105:276-286
  
  # User input
  # Zmean = 14.6 # mean depth
  # Zmax = 25 # max depth
  # Ztherm = 12 # thermocline depth, for graphing purposes
  
  # end user input
  Z = 0:Zmax # range of depths for graphing
  
  p = 6*(Zmean/Zmax)-3
  
  # Depth ratio
  u = Z/Zmax
  u2 = Ztherm/Zmax # for thermocline
  
  # Proportion of lake's surface area that encompasses water dtphs <= Z
  A_u = p*u^2 + (1-p)*u
  A_u2 = p*u2^2 + (1-p)*u2
  
  # Proportion of lakes volume that lies above depth Z
  V_u = (6*u - 3*(1-p)*u^2 - 2*p*u^3) / (3+p)
  V_u2 = (6*u2 - 3*(1-p)*u2^2 - 2*p*u2^3) / (3+p)
  
  # Plots
  if (Plot){
    par(mfrow=c(1,1), mai = c(0.8,0.8, 0.3, 0.1),cex = 0.9)
    plot(Z,A_u,type='l',ylim=c(0,1),xlab = 'Depth (m)', ylab = 'Prop area and volume above Z')
    lines(Z,V_u,col='red')
    abline(v = Ztherm,lty=2)
    abline(h=A_u2,lty=2)
    abline(h=V_u2,lty=2,col='red')
    legend("topleft",legend=c('Area curve','Volume curve','T-cline'),lty=c(1,1,2),col=c('black','red','black'))
  }
  
  Layers = data.frame(AreaAboveZ = A_u2, VolumeAboveZ = V_u2)
  return(Layers)
}
