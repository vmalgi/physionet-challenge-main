function sigP = sigmoidPrime(z)
    sigP = sigmoid(z).*(1-sigmoid(z));
end