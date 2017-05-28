

burn_in <- 1000 # период сходимости
N <- 10000 #кол-во симуляций
getwd()
setwd('C:/Users/asus-пк/Desktop/Учеба/ВШЭ Магистратура/Levy Processes/MCMC Sim/R_Code')
Y <- scan('DailyData.txt')

T <- length(Y) - 2
length(Y)
#-------------------------------------
mu <- rep(NA, burn_in + N)
psi <-rep(NA, burn_in + N)
omega <- rep(NA, burn_in + N)
kappa <- rep(NA, burn_in + N)
theta <- rep(NA, burn_in + N)
lambda <- rep(NA, burn_in + N)
mu_s <- rep(NA, burn_in + N)
sigma2_s <- rep(NA, burn_in + N)
epsilon_s <- rep(NA, (burn_in + N)*T)
epsilon_v <- rep(NA, (burn_in + N)*T)


epsilon_s <- rep(NA, (burn_in + N)*T)
epsilon_v <- rep(NA, (burn_in + N)*T)

# Создаем матрицу столбец, из эксписон на прошлом шаге
epsilon_v <- matrix(epsilon_v, nrow=burn_in + N, ncol = 1)
epsilon_s <- matrix(epsilon_s, nrow=burn_in + N, ncol = 1)
#----------------------------------------------------------
# Введение параметорв для скрытых состояний

V <- rep(NA,(burn_in+N)*(T+2))    # Volatility
V <- matrix(V,(burn_in+N)*(T+2))
Z <- rep(NA,(burn_in+N)*(T+2))   #Size jump
Z <- matrix(Z,(burn_in+N)*(T+2))
B <- rep(NA,(burn_in+N)*(T+2))   # Jump indicate
B <- matrix(B,(burn_in+N)*(T+2))
P <- rep(NA,(burn_in+N)*(T+2))  # Intensivity Lambda jump
P <- matrix(P,(burn_in+N)*(T+2))
#---------------------------------------------------------------------
# Нужно определить гиперпараметры для вспомогательных характеристик

# Hyperparameters for mu 
mu_star <- rep(NA, burn_in +N)
sig_star <- rep(NA, burn_in +N)
mu_star_a <- rep(NA, burn_in +N)
mu_star_b <- rep(NA, burn_in +N)

# Hyperparameters for psi 
psi_star <- rep(NA, burn_in +N)
p_star <- rep(NA, burn_in +N)

# Hyperparameters for omega 
alpha_star <- rep(NA, burn_in +N)
beta_star <- rep(NA, burn_in +N)
beta_star_a <- rep(NA, burn_in +N)
beta_star_b <- rep(NA, burn_in +N)

# Hyperparameters for kappa 
kappa_star <- rep(NA, burn_in +N)
sigma_kappa_star <- rep(NA, burn_in +N)
kappa_star_a <- rep(NA, burn_in +N)
kappa_star_b <- rep(NA, burn_in +N)

# Hyperparameters for theta 
theta_star <- rep(NA, burn_in +N)
sigma_theta_star <- rep(NA, burn_in +N)
theta_star_a <- rep(NA, burn_in +N)
theta_star_b <- rep(NA, burn_in +N)


# Hyperparameters for lambda
alpha_lam <- rep(NA, burn_in +N)
beta_lam <- rep(NA, burn_in +N)

#Hyperparameters for Z (The jump term)
mu_s_star = rep(NA,(burn_in+N)*(T+2))	
mu_s_star = matrix(mu_s_star,nrow=burn_in + N, ncol = 1)
sigma2_s_star = rep(NA,(burn_in+N)*(T+2))		
sigma2_s_star = matrix(sigma2_s_star,nrow=burn_in + N, ncol = 1)

# Закончилась часть для в которой просто ввели переменные.
#--------------------------------------------------------

# Инициализация априорных гиперпараметров
  
mu0 <- 0
sigma0 <- 1
psi0 <- 0
p0 <- 2
alpha_tilde0 <- 2
beta_tilda0 <- 0.005
kappa0 <- 0
sigma_kappa0 <- 1
theta0 <- 0
sigma_theta0 <- 1
alpha_prime <- 2
beta_prime <- 40
mu_s0 <- 0
s0 <- 1
alpha_s <- 5
beta_s <-0.2
delata_t <- 1/252

#----------------------------------------------
# Алгоритм Метрополиса-Хастинга 

pi <- rep(NA, (burn_in+N)*(T+2))
pi <- matrix(pi, nrow = burn_in+N, ncol = 1)
a <- rep(NA, (burn_in+N)*(T+2))
a <- matrix(pi, nrow = burn_in+N, ncol = 1)
accept <- rep(NA, burn_in+N)
#--------------------------------------
#Параметры для первой итерации
mu[1] <- 0.1
omega[1] <- 0.02
psi[1] <- rnorm(1, psi0, sqrt(omega[1]/p0))

kappa[1] <- 5
theta[1] <- 0.15^2
lambda[1] <- rbeta(1, alpha_prime, beta_prime)
mu_s[1] <- 0
sigma2_s[1] <- 0.01
T
P[1,] <- rep(lambda[1], T+2)
P[1,]

