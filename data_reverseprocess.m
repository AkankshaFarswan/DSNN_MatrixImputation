function [X_antilog] = data_reverseprocess(Dmat)
X_antilog = (10.^Dmat) - 1;