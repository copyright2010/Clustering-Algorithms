clear all
num_dimensions = 2;
N = 100; 
%data = [randn(num_dimensions,N)+1  randn(num_dimensions,N)-1];% randn(num_dimensions,N)-8];
data_all = load('data_test');
data = data_all.data_all';

if(num_dimensions ==2)
    plot(data(1,:),data(2,:),'.')
elseif(num_dimensions ==3)
    plot3(data(1,:), data(2,:),data(3,:),'.')
end
pause

num_clusters = 2;
centroids = rand(num_dimensions,num_clusters)*4;

while(1)
    
    distance = zeros(num_clusters,length(data));
    
    %find distance of each data point from each centroids - nned to find
    %the closest centroid
    for k = 1:num_clusters
        sum_dist_squared = zeros(1,length(data));
        for m =1:num_dimensions
            dim_dist(m,:) = data(m,:)-centroids(m,k);
            sum_dist_squared = dim_dist(m,:).^2 +sum_dist_squared; 
        end 
        distance(k,:) = sqrt(sum_dist_squared);
    end
    
    %find the data points closest to each centroid
    %inds will be an integer indicate which if the centroids each data
    %point is closest to
    [vals  inds] = min(distance);
    
    %update centroids
    col_str = 'rgb'; %used to color code the plots
    close all
    orig_centroids = centroids;
    for k = 1:num_clusters
        closest_inds=find(inds==k);
        centroids(:,k) = sum(data(:,closest_inds)'/length(closest_inds));
        if(num_dimensions ==2)
            plot(data(1,closest_inds),data(2,closest_inds),[col_str(k) '.'])
            hold on
            plot(centroids(1,k),centroids(2,k),['k.'],'MarkerSize',40)
            plot(orig_centroids(1,k),orig_centroids(2,k),'kx','MarkerSize',40)       
        elseif(num_dimensions ==3)
            plot3(data(1,closest_inds),data(2,closest_inds),data(3,closest_inds),[col_str(k) '.'])
            hold on
            plot3(centroids(1,k),centroids(2,k),centroids(3,k),['k.'],'MarkerSize',40)
            plot3(orig_centroids(1,k),orig_centroids(2,k),orig_centroids(3,k),'kx','MarkerSize',40)
        end
    end
    pause
end