rm(list=ls(all=TRUE))

options(warn=-1)


library(robustbase)

file <- "C:/temp/electroviewer/mt4R.csv"

MT <- read.csv(file, header=FALSE) # �����l�A���X�g

MT_compoundlist <- MT[,2] # �����ϐ�(���X�g)
MT_sample <- MT[,1] # �ړI�ϐ�(�����l)

# -----------------
#   ���o�X�g��A
# -----------------
brob1 <- lmrob(y~x+I(x^2), data=data.frame(x=MT_compoundlist,y=MT_sample))
p <- predict(brob1)

# ----------------
#   �O���t�̕`��
# ----------------
#plot(MT[,1],MT[,2], xlim=c(0,max(MT[,1],MT[,2])),ylim=c(0,max(MT[,1],MT[,2])))
#par(new=T);plot(c(0:max(MT[,1],MT[,2])), predict(brob1,newdata=data.frame(x=c(0:max(MT[,1],MT[,2])))), type="l",col="blue", lty="solid", xlim=c(0,max(MT[,1],MT[,2])),ylim=c(0,max(MT[,1],MT[,2])), lwd = 1)

# --------------
#   ���ʂ̏o��
# --------------
filepath <- "C:/temp/electroviewer/mt4py.csv"
write.csv(p, filepath, row.names = FALSE, col.names = FALSE)

