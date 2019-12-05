data {
    int<lower=0> T; 
    int<lower=0> dim;
    int<lower=0> k_regimes;
    vector[T] Y;
    matrix[T, dim] X;
}

parameters {
    vector<lower=0, upper=1>[k_regimes] ps;
    vector[k_regimes] consts;
    vector[dim] phis;
    real<lower=0> sigma;
}

model{
    vector[32] stat_probs;
    vector[32] past_probs;
    vector[32] temp_list;
    vector[32] temp_prob;
    matrix[32, 32] prob_matrix;
    vector[4] prob_vec;
    real prob;
    int count;
    int idx;
    real y;
    real denom;
    real likelihood;

    prob_vec[1] = ps[1]; //0->0
    prob_vec[2] = 1-ps[1]; //0->1
    prob_vec[3] = ps[2]; //1->1
    prob_vec[4] = 1-ps[2]; //1->0

    for (k0 in 1 : 4){
        for (k1 in 1 : 4){
            for (k2 in 1 : 4){
                for (k3 in 1 : 4){
                    for (k4 in 1 : 4){
                        count = 1;
                        idx = 1;

                        prob = prob_vec[k0];

                        if ((k0 == 3) || (k0 == 4)) {
                            count += 16;
                        }
                        if ((k0 == 2) || (k0 == 3)) {
                            idx += 16;
                        }

                        if ((k1 == 3) || (k1 == 4)) {
                            count += 8;
                        }
                        if ((k1 == 2) || (k1 == 3)) {
                            idx += 8;
                        }

                        if ((k2 == 3) || (k2 == 4)) {
                            count += 4;
                        }
                        if ((k2 == 2) || (k2 == 3)) {
                            idx += 4;
                        }

                        if ((k3 == 3) || (k3 == 4)) {
                            count += 2;
                        }
                        if ((k3 == 2) || (k3 == 3)) {
                            idx += 2;
                        }

                        if ((k4 == 3) || (k4 == 4)) {
                            count += 1;
                        }
                        if ((k4 == 2) || (k4 == 3)) {
                            idx += 1;
                        }

                        if ((k0==1) && (k1==1)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==1) && (k1==2)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==2) && (k1==3)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==2) && (k1==4)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==3) && (k1==3)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==3) && (k1==4)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==4) && (k1==1)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k0==4) && (k1==2)) {
                            prob = prob*prob_vec[k1];
                        } else{
                            prob = 0;
                        }

                        if ((k1==1) && (k2==1)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k1==1) && (k2==2)) {
                            prob = prob*prob_vec[k2];
                        } else if ((k1==2) && (k2==3)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k1==2) && (k2==4)) {
                            prob = prob*prob_vec[k2];
                        } else if ((k1==3) && (k2==3)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k1==3) && (k2==4)) {
                            prob = prob*prob_vec[k2];
                        } else if ((k1==4) && (k2==1)) {
                            prob = prob*prob_vec[k1];
                        } else if ((k1==4) && (k2==2)) {
                            prob = prob*prob_vec[k2];
                        } else{
                            prob = 0;
                        }

                        if ((k2==1) && (k3==1)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==1) && (k3==2)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==2) && (k3==3)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==2) && (k3==4)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==3) && (k3==3)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==3) && (k3==4)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==4) && (k3==1)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k2==4) && (k3==2)) {
                            prob = prob*prob_vec[k3];
                        } else{
                            prob = 0;
                        }

                        if ((k3==1) && (k4==1)) {
                            prob = prob*prob_vec[k4];
                        } else if ((k3==1) && (k4==2)) {
                            prob = prob*prob_vec[k2];
                        } else if ((k3==2) && (k4==3)) {
                            prob = prob*prob_vec[k3];
                        } else if ((k3==2) && (k4==4)) {
                            prob = prob*prob_vec[k4];
                        } else if ((k3==3) && (k4==3)) {
                            prob = prob*prob_vec[k4];
                        } else if ((k3==3) && (k4==4)) {
                            prob = prob*prob_vec[k4];
                        } else if ((k3==4) && (k4==1)) {
                            prob = prob*prob_vec[k4];
                        } else if ((k3==4) && (k4==2)) {
                            prob = prob*prob_vec[k4];
                        } else{
                            prob = 0;
                        }

                        prob_matrix[count, idx] = prob;
                    }
                }
            }
        }
    }
    
    past_probs = diagonal(prob_matrix);

    for (t in 1:T){
        count = 1;
        for (k0 in 1 : k_regimes){
            for (k1 in 1 : k_regimes){
                for (k2 in 1 : k_regimes){
                    for (k3 in 1 : k_regimes){
                        for (k4 in 1 : k_regimes){
                            y = consts[k0]+phis[k1]*(X[t,1]-consts[k1])+phis[k2]*(X[t,2]-consts[k2])+phis[k3]*(X[t,3]-consts[k3])+phis[k4]*(X[t,4]-consts[k4]);
                            y = exp(-square(Y[t]-y)/(2*sigma));
                            temp_list[count] = y;
                            count += 1;
                        }
                    }
                }
            }
        }

        for (k in 1 : 32){
            if (k == 1){
                denom = temp_list[k];
                likelihood = temp_list[k]*past_probs[k];
            } else {
                denom ＋= temp_list[k];
                likelihood += temp_list[k]*past_probs[k];
            }
        }

        for (k in 1 : 32){
            temp_list[k] = temp_list[k]/denom;
        }

        past_probs = prob_matrix*temp_list;

        target += log(likelihood);
    }
}