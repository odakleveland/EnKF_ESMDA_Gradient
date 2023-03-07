function [Covariance_matrix] = covariance_matrix_m(I,Num_parameters,nr,Pvari,Svari,Rvari)

Covariance_matrix=eye(Num_parameters,Num_parameters); 
for ii=1:length(nr)
    Covariance_matrix(1+I*(ii-1),1+I*(ii-1))=Pvari;
    Covariance_matrix(2+I*(ii-1),2+I*(ii-1))=Svari;
    Covariance_matrix(3+I*(ii-1),3+I*(ii-1))=Rvari;
    
    Covariance_matrix(2+I*(ii-1),1+I*(ii-1))=0; Covariance_matrix(1+I*(ii-1),2+I*(ii-1))=0;
    Covariance_matrix(3+I*(ii-1),1+I*(ii-1))=0;  Covariance_matrix(1+I*(ii-1),3+I*(ii-1))=0; 
    Covariance_matrix(3+I*(ii-1),2+I*(ii-1))=0;  Covariance_matrix(2+I*(ii-1),3+I*(ii-1))=0; 
end

end