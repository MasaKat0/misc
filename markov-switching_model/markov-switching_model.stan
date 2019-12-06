data {
    int<lower=0> T; 
    int<lower=0> k_regimes;
    vector[T] Y;
}

parameters {
    vector<lower=0, upper=1>[k_regimes] ps;
    vector<lower=-2, upper=2>[k_regimes] consts;
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

    prob_matrix[1, 1] = ps[1]; //0->0
    prob_matrix[1, 2] = 1-ps[2]; //1->0
    prob_matrix[2, 1] = 1-ps[1];  //0->1
    prob_matrix[2, 2] = ps[2];  //1->1
    
    past_probs[1] = (1-ps[2])/(2-ps[1]-ps[2]);
    past_probs[2] = (1-ps[1])/(2-ps[1]-ps[2]);

    for (t in 1:T){
        for (k in 1 : k_regimes){
            y = consts[k];
            y = (exp(-square(Y[t]-y)/(2*sigma[k]))/sqrt(2*pi()*sigma[k]))*past_probs[k];
            temp_list[k] = y;
                        
            if (k == 1){
                denom = temp_list[k];
                likelihood = temp_list[k];
            } else {
                denom ＋= temp_list[k];
                likelihood += temp_list[k];
            }
        }

        for (k in 1 : k_regimes){
            temp_list[k] = temp_list[k]/denom;
        }

        past_probs = prob_matrix*temp_list;

        target += log(likelihood);
    }
}