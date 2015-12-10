function [  ] = final( )
    tic
    training = ['Sampledata_1.mat';'Sampledata_2.mat';'Sampledata_3.mat';'Sampledata_4.mat'];
    training = cellstr(training);
    testing = ['Sampledata_test_1.mat';'Sampledata_test_2.mat';'Sampledata_test_3.mat';'Sampledata_test_4.mat'];
    testing = cellstr(testing);
    used_data = testing; %change the variable to training or testing
    files = size(used_data,1);
    for f = 1:files
        tic
        %load Sampledata_1.mat;
        load(used_data{f});
        disp(used_data{f});
        minpts=10;
        a=derivative(spikes);
        comp1 = spike_sort(a,'wav',10);
        fcomp=normalize_vectors(comp1);%normalize_vectors(comp1) 
        mydist=sort(squareform(pdist(fcomp)),1);
        meandist=mean(sort(mydist(minpts+5,:)));
        [out]=dbscan(fcomp,meandist,minpts);
        max(out);
        knn = knn_spikes(fcomp,out,20);

        l = size(out ,1);
        j = 1;
        for i = 1:l
            if(out(i) == -1)
                out(i) = knn(j);
                j = j + 1;
            end
        end
        unique_clusters = unique(out);
        rows = ceil(sqrt(size(unique_clusters,1)+1));
        figure,plot_clustered_spikes(spikes,out);
        plot_clusters(fcomp,out);
        print_stats(out)
        xlswrite(strcat('Results_',num2str(f),'.xlsx'),out);
        toc
    end
end

function [] = print_stats(out)
    unique_clusters = (unique(out));
    cluster_size = size(unique_clusters,1);
    disp(strcat('Cluster size : ',num2str(cluster_size)));
     for i = unique_clusters(1):unique_clusters(cluster_size)
         if i == 0
             continue;
         end
         count = histc(out,i);
         disp(strcat('Cluster count  ',num2str(i),' : ',num2str(count)));
     end
end

function out1 = plot_clusters(comp,out) 
    a=['r' , 'g' ,'b' ,'k' ,'y'];
    unique_clusters = unique(out);
    l = size(unique_clusters,1);
    rows = ceil(sqrt(l+1));
    subplot(rows,rows,rows*rows);
    %sprintf(' for plotiing the points we have shown only first 3 component out of 10')
    %figure('Name','Component Distribution graph')
    hold on 
    title('Component Distribution graph','FontSize',30);
    for i = 1:unique_clusters(l)
        plot3(comp(out==i,1),comp(out==i,2),comp(out==i,3),strcat(a(i),'*'));
    end
    xlabel('wavlet(1)','FontSize',35);ylabel('wavlet(2)','FontSize',35);zlabel('wavlet(3)','FontSize',35);
    hold off
end

function out = plot_clustered_spikes(spikes,out) 
     a=['r' , 'g' ,'b' ,'k' ,'y'];
    unique_clusters = unique(out);
    l = size(unique_clusters,1);
    rows = ceil(sqrt(l+1));
    for i = 1:unique_clusters(l)
     subplot(rows,rows,i);
     %figure('Name',strcat('Spikes of class ',num2str(i)));
     hold on
     plot(1:48,spikes(out==i,:),strcat(a(i),'-'));
     xlabel('Samples','FontSize',35);ylabel('Amplitude','FontSize',35);
     xlim([1,48]);
     title(strcat('Spikes of Cluster : ',num2str(i)),'FontSize',30);
     hold off
    end
end

function [dans]=derivative(spikes)
    [m,n]=size(spikes);
    dans=spikes(:,2:end)-spikes(:,1:end-1);
end

function [nspikes]=normalize_vectors(spikes)
    [m,n]=size(spikes);
    for i=1:n
        tmin=min(spikes(:,i));
        tmax=max(spikes(:,i));
        temp=(spikes(:,i)-tmin)/(tmax-tmin);
        nspikes(:,i)=temp;
    end
end
