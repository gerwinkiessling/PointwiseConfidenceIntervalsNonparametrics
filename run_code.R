###This file produces the tables and figures for my research paper 
#you have to make sure that the working directory of this code is in the folder where the functions.R file is located
rm(list=ls(all=TRUE))
source('functions.R')

# figure 1
plot_density_sim1<-plot_beta_density(2.5, 1.5)
ggsave(filename="plot_density_sim1.png", plot=plot_density_sim1, device="png",width=8, height=6, dpi=500)
#figure2
x_seq=seq(from=0.1, to=0.9, by=0.01)
g_x=g_x_calc(x=x_seq, label='sim1', mu=5/8, var=0.01)
bias<-calculate_bias(label="sim1", x=x_seq , p1=2.5, p2=1.5, bandwidth=0.1, mu=5/8, var=0.01)$bias
E_gh<-g_x+bias
df_plot=data.frame('g_x'=g_x, "E_gh"=E_gh)
plot_bias_sim1<-plot_bias_local_function(df_plot)
ggsave(filename="plot_bias_sim1.png", plot=plot_bias_sim1, device="png",width=8, height=6, dpi=500)

# table 1
undersmooth_seq=c(0, -1/100, -1/50, -1/20, -1/10)
set.seed(1)
#expected value results are in the last column of the coverage matrix
sim1a<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim1", 
                            n=200, p1=2.5, p2=1.5, mu=5/8, var=0.01, x_seq=c(1/6, 5/6), e=-1, r=1,
                            alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="non-robust")
set.seed(2)
sim1b<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim1", 
                           n=400, p1=2.5, p2=1.5, mu=5/8, var=0.01, x_seq=c(1/6, 5/6), e=-1, r=1,
                           alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="non-robust")
sim1b
set.seed(3)
sim1c<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim1", 
                           n=600, p1=2.5, p2=1.5, mu=5/8, var=0.01, x_seq=c(1/6, 5/6), e=-1, r=1,
                           alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="non-robust")
sim1c

# figure 2
x_seq=seq(from=0.1, to=0.9, by=0.01)
g_x=g_x_calc(x=x_seq, label='sim2', mu=5/8, var=0.01)
bias<-calculate_bias(label="sim2", x=x_seq , p1=1.1, p2=1.1, bandwidth=0.1, mu=5/8, var=0.01)$bias
E_gh<-g_x+bias
df_plot=data.frame('g_x'=g_x, "E_gh"=E_gh)
plot_bias_sim2<-plot_bias_linear_function(df_plot)
plot_bias_sim2
ggsave(filename="plot_bias_sim2.png", plot=plot_bias_sim2, device="png",width=8, height=6, dpi=500)

#table 2

set.seed(1)
x_seq<-c(1/4, 3/8, 5/8)
sim2<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim2", 
                     n=200, p1=1.1, p2=1.1, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                     alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="non-robust")
set.seed(2)
x_seq<-c(1/4, 3/8, 5/8)
sim2a<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim2", 
                     n=400, p1=1.1, p2=1.1, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                     alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="non-robust")
#table 3
grid<-matrix(c(1.5, 1, 0.5, 1, 1, 1), 3, 2)

set.seed(1)
replication_multivariate(reps=1000, n=200, factor_undersmooth = -1/100, grid=grid, var=1, cov=0.5, df=4)
set.seed(2)
replication_multivariate(reps=1000, n=200, factor_undersmooth = -1/20, grid=grid, var=1, cov=0.5, df=4)
set.seed(3)
replication_multivariate(reps=1000, n=400, factor_undersmooth = -1/100, grid=grid, var=1, cov=0.5, df=4)
set.seed(4)
replication_multivariate(reps=1000, n=400, factor_undersmooth = -1/20, grid=grid, var=1, cov=0.5, df=4)

#figure 4

x_seq=seq(from=0.01, to=0.99, by=0.0001)
g_x<-g_x_calc(x_seq, "sim4lin", mu=0.01, var=0.01)
df_plot<-data.frame('x_seq'=x_seq, 'g_x'=g_x)
plot_sim4func<-plot_function(df_plot)
ggsave(filename="plot_sim4func.png", plot=plot_sim4func, device="png",width=8, height=6, dpi=500)

#figures 7 and 8


undersmooth_seq=c(0, -1/100, -1/50)
x_seq<-seq(from=0.03, to=0.97, by=0.01)
set.seed(41)
sim4a<-replication_function(reps=1, estimator="local_constant", kernel_function=epa_kernel,label="sim4con", 
                     n=200, p1=1.05, p2=1.05, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                     alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")


set.seed(41)
sim4aa<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim4con", 
                            n=300, p1=1.05, p2=1.05, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                            alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")
sim4aa
set.seed(41)
sim4b<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim4lin", 
                     n=200, p1=1.05, p2=1.05, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                     alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")

set.seed(41)
sim4bb<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim4lin", 
                            n=300, p1=1.05, p2=1.05, mu=5/8, var=0.01, x_seq=x_seq, e=-1, r=1,
                            alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")

plot_mse_function(x_seq=x_seq, data1=sim4a, data2=sim4b, limits=c(0.01, 10), s="sim5_200")
plot_mse_function(x_seq=x_seq, data1=sim4aa, data2=sim4bb, limits=c(0.01, 10), s="sim5_300")
plot_coverage_function(x_seq=x_seq, data1=sim4a, data2=sim42, limits=c(0.2, 1), s="sim4_coverage_200")
plot_coverage_function(x_seq=x_seq, data1=sim4aa, data2=sim4bb, limits=c(0.2, 1), s="sim4_coverage_300")

#figure 5
x_seq=seq(from=0.01, to=0.99, by=0.0001)
g_x<-g_x_calc(x_seq, "sim5lin", mu=0.01, var=0.01)
df_plot<-data.frame('x_seq'=x_seq, 'g_x'=g_x)
plot_sim5func<-plot_function(df_plot)
ggsave(filename="plot_sim5func.png", plot=plot_sim5func, device="png",width=8, height=6, dpi=500)

# figures 9 and 10
undersmooth_seq=c(0, -1/100, -1/50)
x_seq<-seq(from=0.03, to=0.97, by=0.01)
set.seed(5)
sim5a<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim5con", 
                            n=200, p1=1.05, p2=1.05, mu=0.6, var=0.05, x_seq=x_seq, e=-1, r=1,
                            alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")

set.seed(5)
sim5aa<-replication_function(reps=1000, estimator="local_constant", kernel_function=epa_kernel,label="sim5con", 
                             n=300, p1=1.05, p2=1.05, mu=0.6, var=0.05, x_seq=x_seq, e=-1, r=1,
                             alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")
set.seed(5)
sim5b<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim5lin", 
                            n=200, p1=1.05, p2=1.05, mu=0.6, var=0.05, x_seq=x_seq, e=-1, r=1,
                            alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")

set.seed(5)
sim5bb<-replication_function(reps=1000, estimator="local_linear", kernel_function=epa_kernel,label="sim5lin", 
                             n=300, p1=1.05, p2=1.05, mu=0.6, var=0.05, x_seq=x_seq, e=-1, r=1,
                             alpha=0.05, undersmooth_seq=undersmooth_seq, variance_estimator="robust")

plot_mse_function(x_seq=x_seq, data1=sim5a, data2=sim5b, limits=c(0.01, 10), s="sim5_200")
plot_mse_function(x_seq=x_seq, data1=sim5aa, data2=sim5bb, limits=c(0.01, 10), s="sim5_300")
plot_coverage_function(x_seq=x_seq, data1=sim5a, data2=sim5b, limits=c(0.6, 1), s="sim5_coverage_200")
plot_coverage_function(x_seq=x_seq, data1=sim5aa, data2=sim5bb, limits=c(0.6, 1), s="sim5_coverage_300")


####Empirical application#####

###some data cleaning

data<-read.dta("cepr_march_2017.dta")


df<-data.frame(cbind(data['age'],
                     data['female'],
                     data['hrwage'], 
                     data['empl'], 
                     data['educ'] 
))
df<-df[complete.cases(df),]
df<-df[df['educ']=="College", ]
df<-df[df['female']==0, ]
df['experience']=df['age']-22
df=df[df['experience']>0,]
df <- df[df['empl']==1, ]
df <- df[df['age']<=65, ]
df['logwage'] <- log(df['hrwage'])
df <- df[df['logwage']>0,]
set.seed(100)
df<-df[sample(nrow(df), 1000), ]
df

experience=as.numeric(unlist(df['experience']))
logwage=as.numeric(unlist(df['logwage']))
length(experience)
length(logwage)
x_seq=seq(from=min(experience)+2, to=max(experience)-2, length.out=100)

results<-local_linear_empirical(experience, logwage, x_seq)

#figure 6
ggsave(filename="empirical_plot_1.png", plot=results, device="png", dpi=500)
