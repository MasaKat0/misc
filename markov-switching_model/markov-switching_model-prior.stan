data {
    int<lower=0> T; 
    int<lower=0> k_regimes;
    vector[T] Y;
}

parameters {
    vector<lower=0, upper=1>[k_regimes] ps;
    vector[k_regimes] consts;
    vector<lower=0, upper=20>[k_regimes] sigma;
}

model{
    vector[k_regimes] stat_probs;
    vector[k_regimes] past_probs;
    vector[k_regimes] temp_list;
    vector[k_regimes] temp_prob;
    matrix[k_regimes, k_regimes] prob_matrix;
    real y;
    real denom;
    real likelihood;
    int count;
    
    consts[1]~normal(0, 1);
    consts[2]~normal(0, 1);
    sigma[1]~inv_gamma(0.01/2, 0.01/2);
    sigma[2]~inv_gamma(0.01/2, 0.01/2);

    prob_matrix[1, 1] = ps[1]; 
    prob_matrix[1, 2] = 1-ps[1]; 
    prob_matrix[2, 1] = 1-ps[2]; 
    prob_matrix[2, 2] = ps[2]; 
    
    past_probs[1] = ps[1];
    past_probs[2] = 1-ps[1];

    for (t in 1:T){
        count = 1;
        for (k in 1 : k_regimes){
            y = consts[k];
            y = exp(-square(Y[t]-y)/(2*sigma[k]))/sqrt(2*pi()*sigma[k]);
            temp_list[count] = y;
            count += 1;
        }

        for (k in 1 : k_regimes){
            if (k == 1){
                denom = temp_list[k];
                likelihood = temp_list[k]*past_probs[k];
            } else {
                denom ＋= temp_list[k];
                likelihood += temp_list[k]*past_probs[k];
            }
        }

        for (k in 1 : k_regimes){
            temp_list[k] = temp_list[k]/denom;
        }

        past_probs = prob_matrix*temp_list;

        target += log(likelihood);
    }
}