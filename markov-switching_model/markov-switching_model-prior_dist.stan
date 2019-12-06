data {
    int<lower=0> T; 
    int<lower=0> k_regimes;
    vector[T] Y;
}

parameters {
    real<lower=0.00, upper=1.00> prob0;
    real<lower=0.00, upper=1.00> prob1;
    real const0;
    real const1;
    real<lower=0.00> sigma2_0;
    real<lower=0.00> sigma2_1;
}

model{
    vector[k_regimes] past_probs;
    vector[k_regimes] temp_list;

    real denom;
    real likelihood;
    
    const0~normal(0, 1);
    const1~normal(0, 1);
    sigma2_0~inv_gamma(0.01/2, 0.01/2);
    sigma2_1~inv_gamma(0.01/2, 0.01/2);

    past_probs[1] = (1-prob1)/(2-prob0-prob1);
    past_probs[2] = (1-prob0)/(2-prob0-prob1);

    for (t in 1:T){
        temp_list[1] = (exp(-square(Y[t]-const0)/(2*sigma2_0))/sqrt(2*3.14*sigma2_0))*prob0;
        temp_list[2] = (exp(-square(Y[t]-const1)/(2*sigma2_1))/sqrt(2*3.14*sigma2_1))*prob1;
        
        denom = temp_list[1] + temp_list[2];
        likelihood = temp_list[1] + temp_list[2];
        
        temp_list[1] = temp_list[1]/denom;

        past_probs[1] =  temp_list[1]*prob0 + temp_list[2]*(1-prob1);
        past_probs[2] =  temp_list[1]*(1-prob0) + temp_list[2]*prob1;

        target += log(likelihood);
    }
}