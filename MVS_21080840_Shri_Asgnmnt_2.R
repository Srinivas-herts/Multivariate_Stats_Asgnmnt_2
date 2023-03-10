
# Load the required package
library(mlbench)

# Load the Pima Indians Diabetes dataset
data("PimaIndiansDiabetes")

# Create a new dataframe with only the first four columns of the dataset
sample.data <- data.frame(PimaIndiansDiabetes[,1:4])


################ QUESTION 2 a (i): ####################

# Load the required package for Q-Q plotting
library(mvnormtest)

# Create Q-Q plots for each variable
par(mfrow = c(2, 2)) # Set the layout of the plots to 2x2
for (i in 1:4) {
  # Create a Q-Q plot for the ith variable
  qqnorm(sample.data[,i], main = names(sample.data)[i])
  
  # Add a reference line to the Q-Q plot
  qqline(sample.data[,i])
}


################ QUESTION 2 a (ii):  ####################

# Calculate the sample mean vector
mean_vector <- colMeans(sample.data, na.rm = TRUE)
cat("Sample Mean Vector: \n")
print(mean_vector)

# Calculate the sample covariance matrix
cov_matrix <- cov(sample.data, use="pairwise.complete.obs")
cat("\nSample Covariance Matrix: \n")
print(cov_matrix)


################ QUESTION 2 a (iii): ####################

# For multivariate normal distribution functions
library(mvtnorm)
library(ICSNP) # for Hotelling's T-square test

# Define the Hotelling's T-square test statistic
HotellingsT2 <- function(X, mu0, S) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Checking that dimensions of inputs are correct
  if (p != length(mu0)) stop("length of mu0 must equal number of columns of X")
  if (p != dim(S)[1] | p != dim(S)[2]) stop("S must be a square matrix with dimension equal to number of columns of X")
  
  # Calculating the test statistic
  T2 <- n * t(mu0 - colMeans(X)) %*% solve(S) %*% (mu0 - colMeans(X))
  
  return(T2)
}

# Define the null hypothesis mean vector
mu0 <- c(4, 120, 70, 20)

# Calculate the sample mean vector and sample covariance matrix
ybar <- colMeans(sample.data)
S <- cov(sample.data)

# Test the null hypothesis using the Hotelling's T^2 test
T2 <- HotellingsT2(sample.data, mu0, S)
cat("Hotelling's T^2 Test Statistic:", T2, "\n")
cat("Degrees of Freedom:", length(mu0), "\n")
p_value <- pchisq(T2, length(mu0), lower.tail = FALSE)
cat("p-value:", p_value, "\n")

if (p_value < 0.05) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}


################ QUESTION 2 a (iv): ####################

# Plot the sample profile with larger size
plot(sample.data$preg, type = "o", pch = 4, lty = "dashed",
     main = "Pima Indians Diabetes Dataset",
     xlab = "Number of Pregnancies", ylab = "Value", ylim = c(0, 250), cex = 1.5)

lines(sample.data$glu, type = "o", pch = 6, lty = "dotted", col = "red")
lines(sample.data$bp, type = "o", pch = 5, lty = "dotted", col = "green")
lines(sample.data$skin, type = "o", pch = 3, lty = "dotdash", col = "blue")

legend("topright", c("Pregnancies", "Glucose", "Blood Pressure", "Skin Thickness"), 
       lty = c("dashed", "dotted", "dotted", "dotdash"), pch = c(4, 6, 5, 3),
       col = c("black", "red", "green", "blue"), bty = "n")


################ QUESTION 2 a (v): ####################

# install.packages("profileR")
library(profileR)

# Run one sample profile analysis
result1 <- paos(sample.data)

# Print output
print(result1)


## __________ QUESTION 2b _____________



################ QUESTION 2 b (i): #################### 1

sample.data2 <- data.frame(PimaIndiansDiabetes[,c(1:4,9)])

# Split the data into two groups based on diabetes outcome
data_pos <- sample.data2[sample.data2$diabetes == "pos", 1:4]
data_neg <- sample.data2[sample.data2$diabetes == "neg", 1:4]

# Load the required package
library(mvtnorm)

# Calculate sample means and covariance matrices
mean_pos <- apply(data_pos, 2, mean)
mean.neg <- apply(data_neg, 2, mean)
cov_pos <- cov(data_pos)
cov_neg <- cov(data_neg)

# Load the required package
library(Hotelling)

# Perform two-sample T2-test
t2 <- HotellingsT2(rbind(data.pos, data.neg),mu0 = mean.pos - mean.neg, S = (nrow(data.pos) - 1) * cov.pos +
                         (nrow(data.neg) - 1) * cov.neg, n1 = nrow(data.pos), n2 = nrow(data.neg))

t2


################ QUESTION 2 b (ii): #################### 

# to calculate the difference vector
diff.mean <- mean.pos - mean.neg
diff.mean

####################### QUESTION 2 b (iii): ###################### 


# Calculate the pooled covariance matrix
n1 <- nrow(data.pos)
n2 <- nrow(data.neg)
Sp <- ((n1-1)*cov.pos + (n2-1)*cov.neg)/(n1+n2-2)

# Calculate the difference vector between the two sample means
mean.diff <- mean.pos - mean.neg

# Calculate the discriminant function
discriminant <- function(x) {
  x.diff <- x - mean.diff
  return(t(x.diff) %*% solve(Sp) %*% mean.diff)
}

# Apply the discriminant function to the data
discriminant.values <- apply(rbind(data.pos, data.neg), 1, discriminant)

# Compare the discriminant values between the positive and negative groups
summary(discriminant.values)


####################### QUESTION 2 b (iv): ###################### 

# Load the required package
library(profileR)

# Create the profile plot
summary(pbg(sample.data2[,1:4], factor(sample.data2[,5]),
            original.names = TRUE, profile.plot = TRUE))









