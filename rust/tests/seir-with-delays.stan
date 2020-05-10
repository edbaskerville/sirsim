functions {
  real[] ode_ddt(
    real t, real[] state,
    real[] theta,
    real[] x_r, int[] x_i
  ) {
    int n_substates;
    
    int E_gamma_shape;
    int I_gamma_shape;
    int PreD_gamma_shape;
    int OC_gamma_shape;
    int OD_gamma_shape;
    
    real N;
    real b;
    
    real E_mean_duration;
    real I_mean_duration;
    real PreD_mean_duration;
    real OC_mean_duration;
    real OD_mean_duration;
    
    real p_I_R;
    real p_I_PreD;
    
    {
      int index = 1;
      
      n_substates = x_i[index];
      index += 1;
      
      E_gamma_shape = x_i[index];
      index += 1;
      
      I_gamma_shape = x_i[index];
      index += 1;
      
      PreD_gamma_shape = x_i[index];
      index += 1;
      
      OC_gamma_shape = x_i[index];
      index += 1;
      
      OD_gamma_shape = x_i[index];
      index += 1;
    }
    
    {
      int index = 1;
      
      N = x_r[index];
      index += 1;
      
      b = x_r[index];
      index += 1;
      
      E_mean_duration = x_r[index];
      index += 1;
      
      I_mean_duration = x_r[index];
      index += 1;
      
      PreD_mean_duration = x_r[index];
      index += 1;
      
      OC_mean_duration = x_r[index];
      index += 1;
      
      OD_mean_duration = x_r[index];
      index += 1;
      
      p_I_R = x_r[index];
      index += 1;
      
      p_I_PreD = x_r[index];
      index += 1;
    }
    
    {
      real ddt[n_substates];
      
      real S;
      real E[E_gamma_shape];
      real I[I_gamma_shape];
      real R;
      real PreD[PreD_gamma_shape];
      real D;
      real OC[OC_gamma_shape + 1];
      real OD[OD_gamma_shape + 1];
      
      real d_S_E;
      real d_E_I;
      real d_I_R;
      real d_I_PreD;
      real d_PreD_D;
      
      real d_E_E[E_gamma_shape - 1];
      real d_I_I[I_gamma_shape - 1];
      real d_PreD_PreD[PreD_gamma_shape - 1];
      real d_OC_OC[OC_gamma_shape];
      real d_OD_OD[OD_gamma_shape];
      
      real d_E[E_gamma_shape];
      real d_I[I_gamma_shape];
      real d_R;
      real d_PreD[PreD_gamma_shape];
      real d_D;
      real d_OC[OC_gamma_shape + 1];
      real d_OD[OD_gamma_shape + 1];
      
      {
        int index = 1;
        
        E = state[index:(index + E_gamma_shape - 1)];
        index += E_gamma_shape;
        
        I = state[index:(index + I_gamma_shape - 1)];
        index += I_gamma_shape;
        
        R = state[index];
        index += 1;
        
        PreD = state[index:(index + PreD_gamma_shape - 1)];
        index += PreD_gamma_shape;
        
        D = state[index];
        index += 1;
        
        OC = state[index:(index + OC_gamma_shape)];
        index += OC_gamma_shape + 1;
        
        OD = state[index:(index + OD_gamma_shape)];
        index += OD_gamma_shape + 1;
        
        S = N
          - sum(E)
          - sum(I)
          - R
          - sum(PreD)
          - D
        ;
      }
      
      {
        d_S_E = b * (sum(I)) * S / N;
      }
      
      {
        d_E_I = E[E_gamma_shape] * E_gamma_shape / E_mean_duration;
      }
      
      {
        real d_I_ = I[I_gamma_shape] * I_gamma_shape / I_mean_duration;
        d_I_R = p_I_R * d_I_;
        d_I_PreD = p_I_PreD * d_I_;
      }
      
      {
        d_PreD_D = PreD[PreD_gamma_shape] * PreD_gamma_shape / PreD_mean_duration;
      }
      
      for(i in 1:(E_gamma_shape - 1)) {
        d_E_E[i] = E[i] * E_gamma_shape / E_mean_duration;
      }
      
      for(i in 1:(I_gamma_shape - 1)) {
        d_I_I[i] = I[i] * I_gamma_shape / I_mean_duration;
      }
      
      for(i in 1:(PreD_gamma_shape - 1)) {
        d_PreD_PreD[i] = PreD[i] * PreD_gamma_shape / PreD_mean_duration;
      }
      
      for(i in 1:OC_gamma_shape) {
        d_OC_OC[i] = OC[i] * OC_gamma_shape / OC_mean_duration;
      }
      
      for(i in 1:OD_gamma_shape) {
        d_OD_OD[i] = OD[i] * OD_gamma_shape / OD_mean_duration;
      }
      
      if(E_gamma_shape == 1) {
        d_E[1] = d_S_E - d_E_I;
      }
      else {
        d_E[1] = d_S_E - d_E_E[1];
        for(i in 2:(E_gamma_shape - 1)) {
          d_E[i] = d_E_E[i - 1] - d_E_E[i];
        }
        d_E[E_gamma_shape] = d_E_E[E_gamma_shape - 1] - d_E_I;
      }
      
      if(I_gamma_shape == 1) {
        d_I[1] = d_E_I - d_I_R - d_I_PreD;
      }
      else {
        d_I[1] = d_E_I - d_I_I[1];
        for(i in 2:(I_gamma_shape - 1)) {
          d_I[i] = d_I_I[i - 1] - d_I_I[i];
        }
        d_I[I_gamma_shape] = d_I_I[I_gamma_shape - 1] - d_I_R - d_I_PreD;
      }
      
      {
        d_R = d_I_R;
      }
      
      if(PreD_gamma_shape == 1) {
        d_PreD[1] = d_I_PreD - d_PreD_D;
      }
      else {
        d_PreD[1] = d_I_PreD - d_PreD_PreD[1];
        for(i in 2:(PreD_gamma_shape - 1)) {
          d_PreD[i] = d_PreD_PreD[i - 1] - d_PreD_PreD[i];
        }
        d_PreD[PreD_gamma_shape] = d_PreD_PreD[PreD_gamma_shape - 1] - d_PreD_D;
      }
      
      {
        d_D = d_PreD_D;
      }
      
      if(OC_gamma_shape == 0) {
        d_OC[1] = d_E_I;
      }
      else {
        d_OC[1] = d_E_I - d_OC_OC[1];
        for(i in 2:OC_gamma_shape) {
          d_OC[i] = d_OC_OC[i - 1] - d_OC_OC[i];
        }
        d_OC[OC_gamma_shape + 1] = d_OC_OC[OC_gamma_shape];
      }
      
      if(OD_gamma_shape == 0) {
        d_OD[1] = d_PreD_D;
      }
      else {
        d_OD[1] = d_PreD_D - d_OD_OD[1];
        for(i in 2:OD_gamma_shape) {
          d_OD[i] = d_OD_OD[i - 1] - d_OD_OD[i];
        }
        d_OD[OD_gamma_shape + 1] = d_OD_OD[OD_gamma_shape];
      }
      
      {
        int index = 1;
        
        ddt[index:(index + E_gamma_shape - 1)] = d_E;
        index += E_gamma_shape;
        
        ddt[index:(index + I_gamma_shape - 1)] = d_I;
        index += I_gamma_shape;
        
        ddt[index] = d_R;
        index += 1;
        
        ddt[index:(index + PreD_gamma_shape - 1)] = d_PreD;
        index += PreD_gamma_shape;
        
        ddt[index] = d_D;
        index += 1;
        
        ddt[index:(index + OC_gamma_shape)] = d_OC;
        index += OC_gamma_shape + 1;
        
        ddt[index:(index + OD_gamma_shape)] = d_OD;
        index += OD_gamma_shape + 1;
      }
      
      return ddt;
    }
  }
}

data {
  int E_gamma_shape;
  int I_gamma_shape;
  int PreD_gamma_shape;
  int OC_gamma_shape;
  int OD_gamma_shape;
  
  real N;
  real b;
  
  real E_mean_duration;
  real I_mean_duration;
  real PreD_mean_duration;
  real OC_mean_duration;
  real OD_mean_duration;
  
  real p_I_R;
  real p_I_PreD;
  
  real E_init;
  real I_init;
  real R_init;
  real PreD_init;
  real D_init;
  
  int n_times;
  real start_time;
  real times[n_times];
}

transformed data {
  int n_gamma_vars = 5;
  int n_fixed_params = 2;
  int n_probs = 2;
  int n_substates = E_gamma_shape + I_gamma_shape + 1 + PreD_gamma_shape + 1 + OC_gamma_shape + 1 + OD_gamma_shape + 1;
  int x_i[1 + n_gamma_vars];
  real x_r[n_fixed_params + n_gamma_vars + n_probs];
  
  {
    int index = 1;
    
    x_i[index] = n_substates;
    index += 1;
    
    x_i[index] = E_gamma_shape;
    index += 1;
    
    x_i[index] = I_gamma_shape;
    index += 1;
    
    x_i[index] = PreD_gamma_shape;
    index += 1;
    
    x_i[index] = OC_gamma_shape;
    index += 1;
    
    x_i[index] = OD_gamma_shape;
    index += 1;
  }
  
  {
    int index = 1;
    
    x_r[index] = N;
    index += 1;
    
    x_r[index] = b;
    index += 1;
    
    x_r[index] = E_mean_duration;
    index += 1;
    
    x_r[index] = I_mean_duration;
    index += 1;
    
    x_r[index] = PreD_mean_duration;
    index += 1;
    
    x_r[index] = OC_mean_duration;
    index += 1;
    
    x_r[index] = OD_mean_duration;
    index += 1;
    
    x_r[index] = p_I_R;
    index += 1;
    
    x_r[index] = p_I_PreD;
    index += 1;
  }
}

parameters {

}

transformed parameters {

}

model {

}

generated quantities {
  vector[n_times] S;
  vector[n_times] E;
  vector[n_times] I;
  vector[n_times] R;
  vector[n_times] PreD;
  vector[n_times] D;
  vector[n_times] OC;
  vector[n_times] OD;
  
  {
    real initial_state[n_substates];
    real params[0];
    real ode_out[n_times, n_substates];
    
    {
      int index = 1;
      
      initial_state[index:(index + E_gamma_shape - 1)] = rep_array(
        E_init / E_gamma_shape, E_gamma_shape
      );
      index += E_gamma_shape;
      
      initial_state[index:(index + I_gamma_shape - 1)] = rep_array(
        I_init / I_gamma_shape, I_gamma_shape
      );
      index += I_gamma_shape;
      
      initial_state[index] = R_init;
      index += 1;
      
      initial_state[index:(index + PreD_gamma_shape - 1)] = rep_array(
        PreD_init / PreD_gamma_shape, PreD_gamma_shape
      );
      index += PreD_gamma_shape;
      
      initial_state[index] = D_init;
      index += 1;
      
      initial_state[index:(index + OC_gamma_shape)] = rep_array(
        0.0, OC_gamma_shape + 1
      );
      index += OC_gamma_shape + 1;
      
      initial_state[index:(index + OD_gamma_shape)] = rep_array(
        0.0, OD_gamma_shape + 1
      );
      index += OD_gamma_shape + 1;
    }
    
    ode_out = integrate_ode_rk45(
      ode_ddt, initial_state, start_time, times, params, x_r, x_i
    );
    
    for(i in 1:n_times) {
      int index = 1;
      
      E[i] = sum(ode_out[i, index:(index + E_gamma_shape - 1)]);
      index += E_gamma_shape;
      
      I[i] = sum(ode_out[i, index:(index + I_gamma_shape - 1)]);
      index += I_gamma_shape;
      
      R[i] = ode_out[i, index];
      index += 1;
      
      PreD[i] = sum(ode_out[i, index:(index + PreD_gamma_shape - 1)]);
      index += PreD_gamma_shape;
      
      D[i] = ode_out[i, index];
      index += 1;
      
      OC[i] = ode_out[i, index + OC_gamma_shape];
      index += OC_gamma_shape + 1;
      
      OD[i] = ode_out[i, index + OD_gamma_shape];
      index += OD_gamma_shape + 1;
      
      S[i] = N
        - E[i]
        - I[i]
        - R[i]
        - PreD[i]
        - D[i]
      ;
    }
  }
}
