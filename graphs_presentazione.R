PARS <- list(xi=0,alpha=0,omega=1)

H0 <- subset(ALL,delta==PARS$xi & alpha==PARS$alpha & omega==PARS$omega )
righe <- sample(1:nrow(H0),nrow(H0)*.05)

# Build a single long data set with just d and eta1 biases
YY_d <- transform(H0,
                  n = factor(n),
                  ind = "d",
                  bias = Cohens_d - true_d)[, c("n","ind","bias")]

YY_eta <- transform(H0,
                    n = factor(n),
                    ind = "eta",
                    bias = eta1_bias)[, c("n","ind","bias")]

YY_both <- rbind(YY_d, YY_eta)

# One plot: Cohen's d and eta1 biases together
theme_set(theme_bw())
GG_both <- ggplot(YY_both, aes(n, bias, color = ind)) +
  geom_boxplot() +
  theme(legend.title = element_blank()) +
  xlab("sample size") +
  ylab("bias") +
  ggtitle("d = 0     overlap = 1") +
  geom_hline(yintercept = 0, lty = 2)

GG_both


#######

PARS <- list(xi=0,alpha=5,omega=1)

H1 <- subset(ALL,delta==PARS$xi & alpha==PARS$alpha & omega==PARS$omega )
parSN <- pars_skew_normal(alpha=PARS$alpha)

# Build a single long data set with just d and eta1 biases
YY_d <- transform(H1,
                  n = factor(n),
                  ind = "d",
                  bias = Cohens_d - true_d)[, c("n","ind","bias")]

YY_eta <- transform(H1,
                    n = factor(n),
                    ind = "eta",
                    bias = eta1_bias)[, c("n","ind","bias")]

YY_both <- rbind(YY_d, YY_eta)

# One plot: Cohen's d and eta1 biases together
theme_set(theme_bw())
GG_both <- ggplot(YY_both, aes(n, bias, color = ind)) +
  geom_boxplot() +
  theme(legend.title = element_blank()) +
  xlab("sample size") +
  ylab("bias") +
  ggtitle("d = -0.94   overlap = 0.56") +
  geom_hline(yintercept = 0, lty = 2)

GG_both

###########

PARS <- list(xi=0,alpha=0,omega=5)

H2 <- subset(ALL,delta==PARS$xi & alpha==PARS$alpha & omega==PARS$omega )

parSN <- pars_skew_normal(omega=PARS$omega)

# Build a single long data set with just d and eta1 biases
YY_d <- transform(H2,
                  n = factor(n),
                  ind = "d",
                  bias = Cohens_d - true_d)[, c("n","ind","bias")]

YY_eta <- transform(H2,
                    n = factor(n),
                    ind = "eta",
                    bias = eta1_bias)[, c("n","ind","bias")]

YY_both <- rbind(YY_d, YY_eta)

# One plot: Cohen's d and eta1 biases together
theme_set(theme_bw())
GG_both <- ggplot(YY_both, aes(n, bias, color = ind)) +
  geom_boxplot() +
  theme(legend.title = element_blank()) +
  xlab("sample size") +
  ylab("bias") +
  ggtitle("d = 0  overlap = 0.35") +
  geom_hline(yintercept = 0, lty = 2)

GG_both

#######################

library(sn)
library(ggplot2)

# parametri
xi <- 0
omega1 <- 1
omega5 <- 5
alpha0 <- 0
alpha5 <- 5

# griglia
x <- seq(-5, 5, length.out = 500)

# calcolo densità
df <- data.frame(x = x, d1 = d1, d2 = d2)
df_long <- rbind(
  data.frame(x = x, density = d1, dist = "SN(0,1,0)"),
  data.frame(x = x, density = d2, dist = "SN(0,0.5,5)")
)
df$overlap <- pmin(df$d1, df$d2)

# plot
ggplot(df, aes(x, density, color = dist)) +
  geom_line(size = 1) +
  geom_area(data = df, aes(x = x, y = overlap), 
            inherit.aes = FALSE, fill = "green", alpha = 0.3) +
  theme_minimal(base_size = 14) +
  labs(title = "Confronto: Skew-Normal(0,1,0) vs Skew-Normal(0,0,5)",
       x = "y", y = "densità", color = "Distribuzione")


##########

library(sn)
library(ggplot2)

# parameters
xi <- 0
omega1 <- 1
omega5 <- 5
alpha0 <- 0
alpha5 <- 5

# grid
x <- seq(-5, 5, length.out = 1000)

# densities
d1 <- dsn(x, xi = xi, omega = omega1, alpha = alpha0)  # SN(0,1,0)
d2 <- dsn(x, xi = xi, omega = omega1, alpha = alpha5)  # SN(0,omega2,5)

# data frame
df <- data.frame(x = x, d1 = d1, d2 = d2)
df_long <- rbind(
  data.frame(x = x, density = d1, dist = "SN(0,1,0)"),
  data.frame(x = x, density = d2, dist = "SN(0,1,5)")
)

# overlap = minimum of densities
df$overlap <- pmin(df$d1, df$d2)

ggplot(df_long, aes(x, density, color = dist)) +
  geom_line(size = 1) +
  geom_area(data = df, aes(x = x, y = overlap), 
            inherit.aes = FALSE, fill = "green", alpha = 0.3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Confronto: Skew-Normal(0,1,0) vs Skew-Normal(0,1,5)",
    x = "y", y = "densità", color = "Distribuzione"
  ) 

# CASE 2

d1 <- dsn(x, xi = xi, omega = omega1, alpha = alpha0)  # SN(0,1,0)
d2 <- dsn(x, xi = xi, omega = omega1, alpha = alpha5)  # SN(0,omega2,5)

# data frame
df <- data.frame(x = x, d1 = d1, d2 = d2)
df_long <- rbind(
  data.frame(x = x, density = d1, dist = "SN(0,1,0)"),
  data.frame(x = x, density = d2, dist = "SN(0,5,0)")
)

# overlap = minimum of densities
df$overlap <- pmin(df$d1, df$d2)

ggplot(df_long, aes(x, density, color = dist)) +
  geom_line(size = 1) +
  geom_area(data = df, aes(x = x, y = overlap), 
            inherit.aes = FALSE, fill = "green", alpha = 0.3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Confronto: Skew-Normal(0,1,0) vs Skew-Normal(0,1,5)",
    x = "y", y = "densità", color = "Distribuzione"
  ) 


## absolute mean bias

(BIAS_eta1 <- aggregate(eta1_bias~n*delta*alpha*omega,data = ALL,FUN=mean))
colnames(BIAS_eta1) <- gsub("eta1_","",colnames(BIAS_eta1))
BIAS_eta1$V <- aggregate(eta1_bias~n*delta*alpha*omega,data = ALL,FUN=var)$eta1_bias
BIAS_eta1$index <- "eta"

BIAS_cohen <- aggregate(Cohens_d_bias~n*delta*alpha*omega,data = ALL,FUN=mean)
colnames(BIAS_cohen) <- gsub("Cohens_d_","",colnames(BIAS_cohen))
BIAS_cohen$V <- aggregate(Cohens_d_bias~n*delta*alpha*omega,data = ALL,FUN=var)$Cohens_d_bias
BIAS_cohen$index <- "Cohen's d"

BIAS <- rbind(BIAS_eta1,BIAS_cohen)
BIAS$n <- factor(BIAS$n)
BIAS$x <- as.numeric(BIAS$n)
rm(BIAS_eta1,BIAS_cohen)

ggplot(subset(BIAS,delta==0),aes(x,bias,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_hline(yintercept = 0, lty = 2)+geom_point()+geom_line()+xlab("sample size")+ylab("absolute mean bias")+scale_x_continuous(breaks = unique(BIAS$x),labels = levels(BIAS$n))

ggplot(subset(BIAS,delta==1),aes(x,bias,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_hline(yintercept = 0, lty = 2)+geom_point()+geom_line()+xlab("sample size")+ylab("absolute mean bias")+scale_x_continuous(breaks = unique(BIAS$x),labels = levels(BIAS$n))

## relative mean bias

RBIAS_eta1 <- aggregate(eta1_rbias~n*delta*alpha*omega,data = ALL,FUN=mean)
colnames(RBIAS_eta1) <- gsub("eta1_","",colnames(RBIAS_eta1))
RBIAS_eta1$index <- "eta"


RBIAS_cohen <- aggregate(Cohens_d_rbias~n*delta*alpha*omega,data = ALL,FUN=mean)
colnames(RBIAS_cohen) <- gsub("Cohens_d_","",colnames(RBIAS_cohen))
RBIAS_cohen$index <- "Cohen's d"

RBIAS <- rbind(RBIAS_eta1,RBIAS_cohen)
RBIAS$n <- factor(RBIAS$n)
RBIAS$x <- as.numeric(RBIAS$n)
rm(RBIAS_eta1,RBIAS_cohen)

ggplot(subset(RBIAS,delta==0),aes(x,rbias,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_hline(yintercept = c(-.1,.1), lty = 2)+geom_point()+geom_line()+xlab("sample size")+ylab("relative mean bias")+scale_x_continuous(breaks = unique(RBIAS$x),labels = levels(RBIAS$n))

ggplot(subset(RBIAS,delta==1),aes(x,rbias,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_hline(yintercept = c(-.1,.1), lty = 2)+geom_point()+geom_line()+xlab("sample size")+ylab("relative mean bias")+scale_x_continuous(breaks = unique(BIAS$x),labels = levels(BIAS$n))


# Mean square error

ALL$eta1_bias2 <- ALL$eta1_bias^2
ALL$Cohens_d_bias2 <- ALL$Cohens_d_bias^2

BIAS2_eta1 <- aggregate(eta1_bias2~n*delta*alpha*omega,data = ALL,FUN=mean)
colnames(BIAS2_eta1) <- gsub("eta1_","",colnames(BIAS2_eta1))
BIAS2_eta1$V <- aggregate(eta1_bias~n*delta*alpha*omega,data = ALL,FUN=var)$eta1_bias
BIAS2_eta1$index <- "eta"

BIAS2_cohen <- aggregate(Cohens_d_bias2~n*delta*alpha*omega,data = ALL,FUN=mean)
colnames(BIAS2_cohen) <- gsub("Cohens_d_","",colnames(BIAS2_cohen))
BIAS2_cohen$V <- aggregate(Cohens_d_bias~n*delta*alpha*omega,data = ALL,FUN=var)$Cohens_d_bias
BIAS2_cohen$index <- "Cohen's d"

BIAS2 <- rbind(BIAS2_eta1,BIAS2_cohen)
BIAS2$n <- factor(BIAS2$n)
BIAS2$x <- as.numeric(BIAS2$n)
rm(BIAS2_eta1,BIAS2_cohen)
BIAS <- BIAS2
BIAS$MSE <- BIAS$bias2

ggplot(subset(BIAS,delta==0),aes(x,MSE,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_point()+geom_line()+xlab("sample size")+ylab("mean square error")+scale_x_continuous(breaks = unique(RBIAS$x),labels = levels(RBIAS$n))
ggplot(subset(BIAS,delta==1),aes(x,MSE,color=index,shape=index))+theme_bw()+facet_grid(alpha~omega)+geom_point()+geom_line()+xlab("sample size")+ylab("mean square error")+scale_x_continuous(breaks = unique(RBIAS$x),labels = levels(RBIAS$n))
