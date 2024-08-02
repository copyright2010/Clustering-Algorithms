num_dimensions = 2;
N = 1000; 
data = [randn(num_dimensions,N)+2  randn(num_dimensions,N)-2];
% if(num_dimensions == 2)
%     plot(data(1,:),data(2,:),'.')
% elseif(num_dimensions ==3)
%     plot3(data(1,:), data(2,:),data(3,:),'.')
% end

% data = data';
% 
% cov_matrix = [0.1 0.8 ; 0.5 1.27];
% 
% x1 = data(1:N,1)';
% x2 = data(1:N,2)';
% x3 = x1+12;
% x4 = x2-5;
% 
% data_new = cov_matrix*[x1;x2];
% data_new2 = cov_matrix*[x3;x4];
% data_all = [data_new data_new2];
% figure(2)
% plot(data_all(1,:),data_all(2,:), '.')
% data_all = data_all';

data_all = load('data_test');
data_all = data_all.data_all;
figure(2)
plot(data_all(1,:),data_all(2,:), '.')
C = 0;
epsilon = 0.4;
MinPts = 45;
n = size(data_all,1);
IDX = zeros(n,1);

D = pdist2(data_all,data_all);

visited=false(n,1); % set visited as an array of 0's - same size of data
isnoise=false(n,1); %

for i=1:n
    if ~visited(i) % if increment of visited is 0
        visited(i)=true; % 
        
        Neighbours = find(D(i,:)<=epsilon);
        
        if numel(Neighbours)<MinPts
            isnoise(i)=true;
        else
            C=C+1;
            IDX(i) = C;
            
            k = 1;
            while true
                j = Neighbours(k);
                
                if ~visited(j)
                    visited(j)=true;
                    Neighbours2 = find(D(j,:)<=epsilon);
                    if numel(Neighbours2)>=MinPts
                        Neighbours = [Neighbours Neighbours2];
                    end
                end
                if IDX(j) == 0
                    IDX(j) = C;
                end
                
                k = k + 1;
                if k > numel(Neighbours)
                    break;
                end
            end
        end
    end
end

PlotClusterinResult(data_all,IDX)
function PlotClusterinResult(X, IDX)

    k=max(IDX);

    Colors=hsv(k);

    Legends = {};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = 'x';
            MarkerSize = 8;
            Color = Colors(i,:);
            Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 6;
            Color = [0 0 0];
            if ~isempty(Xi)
                Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            plot(Xi(:,1),Xi(:,2),'.')%,Style,'MarkerSize',MarkerSize,'Color',Color);
        end
        hold on;
    end
    hold off;
    %axis equal;
    %grid on;
    %legend(Legends);
    %legend('Location', 'NorthEastOutside');

end


































% 
% dist = 0.5;
% 
% x = data(1,:);
% y = data(2,:);
% 
% j = 2;
% for i = 1:length(x)-1
%     pdist_x(i) = pdist2(x(i),x(j));
%     j = j + 1;
% end
% 
% j = 2;
% for i = 1:length(y)-1
%     pdist_y(i) = pdist2(y(i),y(j));
%     j = j + 1;
% end
% 
% n = length(x)*length(y);
% 
% for k = 1:length(x)
%     for i = 1:length(x)
%         pdist = [pdist; pdist2(x(k),x(i))];
%     end
% end