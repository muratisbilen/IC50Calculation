setwd("inputFiles")
inputFiles=try(shell("dir /B ",intern=T,wait=T))
setwd("../")

sigmoid3 = function(params, x) {
	params[1]+((params[2]-params[1])/(1+10^(x-params[3])))
}
sigmoid3B = function(params, x) {
	0+((params[1]-0)/(1+10^(x-params[2])))
}
sigmoid3T = function(params, x) {
	params[1]+((100-params[1])/(1+10^(x-params[2])))
}

sigmoid4 = function(params, x) {
	params[1]+((params[2]-params[1])/(1+10^((params[3]-x)*params[4])))
}
sigmoid4B = function(params, x) {
	0+((params[1]-0)/(1+10^((params[2]-x)*params[3])))
}
sigmoid4T = function(params, x) {
	params[1]+((100-params[1])/(1+10^((params[2]-x)*params[3])))
}

mat = c("Model","Sample","Drug","log(IC50)","Std. Err of Model","Std. Err of log(EC50)",
	"log(EC50)","Activity Area","Amax","log(IC90)","log(IC95)")

for(files in inputFiles){
	outputFile = paste(substr(files,1,nchar(files)-4)," Drug Data.txt",sep="")
	outputPDF = paste(substr(files,1,nchar(files)-4)," Fitted Curves.pdf",sep="")
	pdf(outputPDF,width=14)
	setwd("inputFiles")
	d=read.delim(files, sep="\t",header=F)
	setwd("../")
	pb <- winProgressBar(title=substr(files,1,nchar(files)-4), label="0% done", min=0, max=100, initial=0)
	x1=log10(as.numeric(as.matrix(d[1,-c(1,2)])))
	x2 = seq(min(x1,na.rm=T),max(x1,na.rm=T),0.01)
	for(j in 2:nrow(d)){
		drug = as.vector(as.matrix(d[j,2]))
		cell = as.vector(as.matrix(d[j,1]))
		l3="3-Parameter"
		l3B="3-Parameter Bottom 0"
		l3T="3-Parameter Top 100"
		l4="4-Parameter"
		l4B="4-Parameter Bottom 0"
		l4T="4-Parameter Top 100"

		y=as.numeric(as.matrix(d[j,-c(1,2)]))
		if(is.vector(y)){
			y = cbind(y,y)
		}
		ym1 = apply(y,1,mean,na.rm=T)
		
		l = c()
		for(i in 1:nrow(y)){
			if(length(y[i,which(!is.na(y[i,]))])>1){
				l = c(l,sd(y[i,],na.rm=T)/sqrt(length(which(!is.na(y[i,])))))
			}else{
				l = c(l,0)
			}
		}
		
		ylmax = max(ym1+l,na.rm=T)
		ylmin = min(ym1-l,na.rm=T)
		
		par(mfrow = c(2,3))
		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit3 = nls(ym~Ab+(At-Ab)/(1+10^(x-EC50)), start=list(Ab = min(ym), At=max(ym), EC50=0))
			params3 = coef(fit3)
			ic3=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(50-coef(fit3)[1])-1)
			ic3.90=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(10-coef(fit3)[1])-1)
			ic3.95=coef(fit3)[3]+log10((coef(fit3)[2]-coef(fit3)[1])/(5-coef(fit3)[1])-1)
			icse3 = summary(fit3)$sigma
			yp3 = sigmoid3(params3,x2)
			aa3 = sum(0.01*(params3[2]-yp3))
			amax3 = params3[2]-params3[1]
			ecse3 = summary(fit3)$coefficients[3,2]
			
			
			plot(x2,yp3,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3,0,na.rm=T),max(ylmax,yp3,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3,"\nError: ",round(icse3,4),"\tlog(IC50): ",round(ic3,4),sep=""))
			
			
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic3,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic3<<-NA;ic3.90<<-NA;ic3.95<<-NA;icse3<<-NA;aa3<<-NA;amax3<<-NA;yp3<<-0;
			params3 <<- c(NA,NA,NA); ecse3 <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,0,na.rm=T),max(ylmax,yp3,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3,sep=""))

			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})
		
		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit3B = nls(ym~0+(At-0)/(1+10^(x-EC50)), start=list(At=max(ym), EC50=0))
			params3B = coef(fit3B)
			ic3B=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(50-0)-1)
			ic3B.90=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(10-0)-1)
			ic3B.95=coef(fit3B)[2]+log10((coef(fit3B)[1]-0)/(5-0)-1)
			icse3B = summary(fit3B)$sigma
			yp3B = sigmoid3B(params3B,x2)
			aa3B = sum(0.01*(params3B[1]-yp3B))
			amax3B = params3B[1]-0
			ecse3B = summary(fit3B)$coefficients[2,2]
			
			plot(x2,yp3B,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3B,0,na.rm=T),max(ylmax,yp3B,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3B,"\nError: ",round(icse3B,4),"\tlog(IC50): ",round(ic3B,4),sep=""))
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic3B,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic3B<<-NA;ic3B.90<<-NA;ic3B.95<<-NA;icse3B<<-NA;aa3B<<-NA;amax3B<<-NA;yp3B<<-0;
			params3B <<- c(NA,NA);ecse3B <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp3B,0,na.rm=T),max(ylmax,yp3B,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3B,sep=""))
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})

		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit3T = nls(ym~Ab+(100-Ab)/(1+10^(x-EC50)), start=list(Ab=min(ym), EC50=0))
			params3T = coef(fit3T)
			ic3T=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(50-coef(fit3T)[1])-1)
			ic3T.90=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(10-coef(fit3T)[1])-1)
			ic3T.95=coef(fit3T)[2]+log10((100-coef(fit3T)[1])/(5-coef(fit3T)[1])-1)
			icse3T = summary(fit3T)$sigma
			yp3T = sigmoid3T(params3T,x2)
			aa3T = sum(0.01*(100-yp3T))
			amax3T = 100-params3T[1]
			ecse3T = summary(fit3T)$coefficients[2,2]
			
			plot(x2,yp3T,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp3T,0,na.rm=T),max(ylmax,yp3T,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3T,"\nError: ",round(icse3T,4),"\tlog(IC50): ",round(ic3T,4),sep=""))
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic3T,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic3T<<-NA;ic3T.90<<-NA;ic3T.95<<-NA;icse3T<<-NA;aa3T<<-NA;amax3T<<-NA;yp3T<<-0;
			params3T <<- c(NA,NA);ecse3T <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp3T,0,na.rm=T),max(ylmax,yp3T,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l3T,sep=""))
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})
		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit4 = nls(ym~Ab+(At-Ab)/(1+10^((EC50-x)*H)), start=list(Ab = min(ym), At=max(ym), EC50=0,H=-1))
			params4 = coef(fit4)
			ic4=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(50-coef(fit4)[1])-1)/coef(fit4)[4]
			ic4.90=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(10-coef(fit4)[1])-1)/coef(fit4)[4]
			ic4.95=coef(fit4)[3]-log10((coef(fit4)[2]-coef(fit4)[1])/(5-coef(fit4)[1])-1)/coef(fit4)[4]
			icse4 = summary(fit4)$sigma
			yp4 = sigmoid4(params4,x2)
			aa4 = sum(0.01*(params4[2]-yp4))
			amax4 = params4[2]-params4[1]
			ecse4 = summary(fit4)$coefficients[3,2]
			
			plot(x2,yp4,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4,0,na.rm=T),max(ylmax,yp4,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4,"\nError: ",round(icse4,4),"\tlog(IC50): ",round(ic4,4),sep=""))
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic4,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic4<<-NA;ic4.90<<-NA;ic4.95<<-NA;icse4<<-NA;aa4<<-NA;amax4<<-NA;yp4<<-0;
			params4 <<- c(NA,NA,NA,NA); ecse4 <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4,0,na.rm=T),max(ylmax,yp4,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4,sep=""))
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})
		
		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit4B = nls(ym~0+(At-0)/(1+10^((EC50-x)*H)), start=list(At=max(ym), EC50=0,H=-1))
			params4B = coef(fit4B)
			ic4B=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(50-0)-1)/coef(fit4B)[3]
			ic4B.90=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(10-0)-1)/coef(fit4B)[3]
			ic4B.95=coef(fit4B)[2]-log10((coef(fit4B)[1]-0)/(5-0)-1)/coef(fit4B)[3]
			icse4B = summary(fit4B)$sigma
			yp4B = sigmoid4B(params4B,x2)
			aa4B = sum(0.01*(params4B[1]-yp4B))
			amax4B = params4B[1]-0
			ecse4B = summary(fit4B)$coefficients[2,2]
			
			plot(x2,yp4B,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4B,0,na.rm=T),max(ylmax,yp4B,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4B,"\nError: ",round(icse4B,4),"\tlog(IC50): ",round(ic4B,4),sep=""))
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic4B,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic4B<<-NA;ic4B.90<<-NA;ic4B.95<<-NA;icse4B<<-NA;aa4B<<-NA;amax4B<<-NA;yp4B<<-0;
			params4B <<- c(NA,NA,NA); ecse4B <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4B,0,na.rm=T),max(ylmax,yp4B,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4B,sep=""))
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})
		
		tryCatch({
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			fit4T = nls(ym~Ab+(100-Ab)/(1+10^((EC50-x)*H)), start=list(Ab = min(ym), EC50=0,H=-1))
			params4T = coef(fit4T)
			ic4T=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(50-coef(fit4T)[1])-1)/coef(fit4T)[3]
			ic4T.90=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(10-coef(fit4T)[1])-1)/coef(fit4T)[3]
			ic4T.95=coef(fit4T)[2]-log10((100-coef(fit4T)[1])/(5-coef(fit4T)[1])-1)/coef(fit4T)[3]
			icse4T = summary(fit4T)$sigma
			yp4T = sigmoid4T(params4T,x2)
			aa4T = sum(0.01*(100-yp4T))
			amax4T = 100-params4T[1]
			ecse4T = summary(fit4T)$coefficients[2,2]
			
			plot(x2,yp4T,type="l",lwd=4,col="red",ylim=c(min(ylmin,yp4T,0,na.rm=T),max(ylmax,yp4T,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4T,"\nError: ",round(icse4T,4),"\tlog(IC50): ",round(ic4T,4),sep=""))
			points(x,ym,pch=16,cex=1.5)
			abline(h=50,v=ic4T,lty=2)
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		}, error=function(e){
			posym = which(!is.na(ym1))
			ym = ym1[posym]
			x = x1[posym]
			ic4T<<-NA;ic4T.90<<-NA;ic4T.95<<-NA;icse4T<<-NA;aa4T<<-NA;amax4T<<-NA;yp4T<<-0;
			params4T <<- c(NA,NA,NA); ecse4T <<- NA
			plot(x,ym,pch=16,cex=1.5,ylim=c(min(ylmin,yp4T,0,na.rm=T),max(ylmax,yp4T,na.rm=T)),
				xlab=expression(paste("Concentration (log(",mu,"M))")), ylab="Growth (%)",
				main=paste(cell," - ",drug,"\n",l4T,sep=""))
			for(i in 1:nrow(y)){
				lines(c(x[i],x[i]),c(ym[i]-l[i],ym[i]+l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]-l[i],ym[i]-l[i]),lwd=2)
				lines(c(x[i]-0.05,x[i]+0.05),c(ym[i]+l[i],ym[i]+l[i]),lwd=2)
			}
		})
		
		stanerr = c(icse3,icse3B,icse3T,icse4,icse4B,icse4T)
		minpos = which(stanerr == min(stanerr,na.rm=T))
		adjs = c(0.03, 0.367, 0.705, 0.03, 0.367, 0.705)
		padjs = c(2, 2, 2, 11.7, 11.7, 11.7)
		
		if(length(minpos)>0){
			mtext("*", adj=adjs[minpos],padj=padjs[minpos],outer=T,col="blue",cex=3)
			mtext("*", adj=adjs[minpos]+0.285,padj=padjs[minpos],outer=T,col="blue",cex=3)
		}	
		
		r_3 = c(l3,cell,drug,ic3,icse3,ecse3,params3[3],aa3,amax3,
			ic3.90,ic3.95)
		r_3B = c(l3B,cell,drug,ic3B,icse3B,ecse3B,params3B[2],aa3B,
			amax3B,ic3B.90,ic3B.95)
		r_3T = c(l3T,cell,drug,ic3T,icse3T,ecse3T,params3T[2],aa3T,
			amax3T,ic3T.90,ic3T.95)
		r_4 = c(l4,cell,drug,ic4,icse4,ecse4,params4[3],aa4,amax4,
			ic4.90,ic4.95)
		r_4B = c(l4B,cell,drug,ic4B,icse4B,ecse4B,params4B[2],aa4B,
			amax4B,ic4B.90,ic4B.95)
		r_4T = c(l4T,cell,drug,ic4T,icse4T,ecse4T,params4T[2],aa4T,
			amax4T,ic4T.90,ic4T.95)
		
		
		mat = rbind(mat,r_3,r_3B,r_3T,r_4,r_4B,r_4T)
		
		info <- sprintf("%d%% done", round((j/nrow(d))*100))
		setWinProgressBar(pb, (j/nrow(d))*100, label=info)
	}
	close(pb)
	dev.off()
	write(file=outputFile,t(mat),ncol=ncol(mat),sep="\t")
}

