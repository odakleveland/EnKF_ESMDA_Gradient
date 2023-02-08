function obs1 = reflectivity_convolution(time,theta,Rpp,w)
%Convolution, and observed data
obs1=zeros(1,((length(time)-1)*length(theta)));
for kk=1:length(theta)
    con=conv(Rpp(:,kk),w, 'same'); %USE METHOD SAME ???
    obs1((kk-1)*length(time)+1:kk*length(time))=con;
end
end