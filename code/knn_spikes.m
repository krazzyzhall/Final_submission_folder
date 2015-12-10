function [ output ] = knn_spikes( component, input, nn)
    unique_clusters = unique(input);
    cluster_size = size(unique_clusters,1);
    for i=1:cluster_size -1
        c{i}=component(input==i,:);
        [l,m]=size(c{i});
        dist{i} = zeros(l,1);
        counter{i} = 0;
        index{i} = 1;
    end
    temp = component(input==-1,:);
    [l1,m1] = size(temp);
    output = zeros(l1,1);
    for i = 1:l1
        for j = 1:cluster_size -1
            a = pdist([temp(i,:) ; c{j}], 'euclidean');
            a =  a(1:size(c{j},1));
            %dist{j} = pdist([repmat(temp(i,:),size(c{j},1),1) ; c{j}], 'euclidean');
            dist{j} = sort(a);
        end
        while(add_cell(counter,cluster_size) < nn)
            minimum = min_cell(dist , index,cluster_size);
            index{minimum} = index{minimum} + 1;
            counter{minimum} = counter{minimum} + 1;
        end
        output(i) = max_cell(counter,cluster_size);
        for j = 1:cluster_size - 1;
            counter{j} = 0;
            index{j} = 1;
        end
    end
end


function [out] = add_cell(C , cluster_size)
    out = 0;
    for i = 1:cluster_size - 1
        out = out + C{i};
    end
end

function [out] = min_cell(C , index,cluster_size)
    min = 9999;
    for i = 1:cluster_size - 1
        if(size(C{i},1) > index{i})
            continue;
        end
        if(C{i}(index{i}) < min)
            min = C{i}(index{i});
            out = i;
        end
    end
end

function [out] = max_cell(counter,cluster_size)
    max = 0;
    for i = 1:cluster_size - 1
        if(counter{i} > max)
            max = counter{i};
            out = i;
        end
    end
end
