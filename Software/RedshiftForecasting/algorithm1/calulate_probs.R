#########
######### by James Long
######### date Nov 26, 2010
#########
######### calculate the probability of a feature other than uvot detection
######### being selected by cart
#########


edge_prob = function(n11,n12,n21,n22){
  first = factorial(n11 + n21) / (factorial(n11) * factorial(n21))
  second = factorial(n12 + n22) / (factorial(n12) * factorial(n22))
  third = (factorial(n11 + n12) * factorial(n21 + n22)) / factorial(n11 + n12 + n21 + n22)
  answer = 2 * first * second * third
  return(answer)
}


n11 = 11
n12 = 3
n21 = 31
n22 = 81
num_vars = 40
the_prob = edge_prob(n11,n12,n21,n22) + edge_prob(n11+1,n12-1,n21,n22) + edge_prob(n11+2,n12-2,n21,n22) + edge_prob(n11+3,n12-3,n21,n22)
#the_prob = edge_prob(n11,n12,n21,n22) + edge_prob(n11+1,n12-1,n21,n22) + edge_prob(n11+2,n12-2,n21,n22) + edge_prob(n11+3,n12-3,n21,n22) + edge_prob(n11+3,n12-3,n21,n22)


(1 - the_prob) ^ num_vars
