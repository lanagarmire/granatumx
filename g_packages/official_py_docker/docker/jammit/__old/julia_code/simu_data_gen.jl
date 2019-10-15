function simu_dat_gen(
  ;
  n_features=20,
  n_samples=10,
  n_sig_features=5,
  sig_intensity=10,
  noise_intensity=1,
  features_as_cols=true,
)
  if n_sig_features > n_features
    error("n_sig_features ($n_sig_features) > n_features ($n_features)")
  end

  mat = randn(n_features, n_samples)

  if features_as_cols
    return mat
  else
    return mat'
  end
end


# test

print("Hello!")