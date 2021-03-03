%% Loading Data and initialization
clear all
close all
clc

%% Adding every files in the path

addpath(genpath(pwd))


datasets = {"indian"}; 
options.lambda=0.001;
options.delta=0.01;
options.display = true;
options.inneriter = 10;
options.delta_iter = 1e-2;
options.save = false;
n_samples = 300; % If options.save == true, number of points in the .dat file
n_runs = 20; % Number of runs tested per dataset

methods = {"titan\_minvol","minvol"};
symbols = {'-o','-s'};

for dataset = datasets
	%% dataset loading
    
    [X,r,maxtime] = dataset_loader(dataset{:});
    options.maxtime = maxtime;
	k = 0;
	nX = norm(X,'fro')^2;
	[m,n] = size(X);
    methods_computed = {};
	for run = 1:n_runs
        %% Initialization
        
		W = rand(m,r);
		H = FGMfcnls(X,W,[],100);   

		scale = sum(W);
		W = W./scale;
		H = H.*scale';

		options.W = W;
		options.H = H;

		%% Running the methods

		% TITANized-minvol
		if ~isempty(find([methods{:}] == "titan\_minvol"))
			k = k+1;
			disp("computing TITANized-minvol... "+dataset{:}+"... run "+num2str(run))
			options.inertial = true;
			[W,H,e,en,ep,etx,lambda] = titanminvol(X,r,options);
			es{k} = e;
			ens{k} = en;
			eps{k} = ep;
			etxs{k} = etx;
            if run==1
                methods_computed{k} = "titan\_minvol";
            end
		end

		% halsiminvolNMF
		% Without extrapolation
		if ~isempty(find([methods{:}] == "hals"))
			k = k+1;
			disp("computing hals... "+dataset{:}+"... run "+num2str(run))
			options.inertial = false;
			[W,H,e,en,ep,etx,lambda] = halsiminvolNMF(X,r,options);
			es{k} = e;
			ens{k} = en;
			eps{k} = ep;
			etxs{k} = etx;
            if run==1
                methods_computed{k} = "hals";
            end
		end

		% With extrapolation
		if ~isempty(find([methods{:}] == "hals\_extra"))
			k = k+1;
			disp("computing hals\_extra... "+dataset{:}+"... run "+num2str(run))
			options.inertial = true;
			[W,H,e,en,ep,etx,lambda] = halsiminvolNMF(X,r,options); 
			es{k} = e;
			ens{k} = en;
			eps{k} = ep;
			etxs{k} = etx;
            if run==1
                methods_computed{k} = "hals\_extra";
            end
		end

		% minvolNMF
		if ~isempty(find([methods{:}] == "minvol"))
			k = k+1;
			disp("computing minvol... "+dataset{:}+"... run "+num2str(run))
			[W,H,e,en,ep,etx,lambda] = minvolNMF(X,r,options);
			es{k} = e;
			ens{k} = en;
			eps{k} = ep;
			etxs{k} = etx;
            if run==1
                methods_computed{k} = "minvol";
            end
		end
	end

	%% Homogenizing time sampling
	es_min = min([es{:}]);
	ens_min = min([ens{:}]);
	eps_min = min([eps{:}]);
	minmaxtime = options.maxtime;
	for kk = 1:k
		minmaxtime = min(minmaxtime,max(etxs{kk}));
	end
	new_time_sampling = 0:(minmaxtime/(n_samples-1)):minmaxtime;

	for kk = 1:k
		sampled_es{kk} = interp1(etxs{kk},(es{kk}-es_min)./nX,new_time_sampling)';
		sampled_ens{kk} = interp1(etxs{kk},(ens{kk}-ens_min)./nX,new_time_sampling)';
		sampled_eps{kk} = interp1(etxs{kk},(eps{kk}-eps_min)./nX,new_time_sampling)';
	end
	
	%% Computing the average for the loaded dataset and saving results
	Means = {};
	for kk = 1:length(methods)
		Means{1,kk} = mean([sampled_es{kk:length(methods):k}],2);
		Means{2,kk} = mean([sampled_ens{kk:length(methods):k}],2);
		Means{3,kk} = mean([sampled_eps{kk:length(methods):k}],2);
		if options.save
			todat = cat(2,new_time_sampling',Means{1,kk});
			save(pwd+"\saved_data\mean_"+dataset+num2str(kk)+".dat", 'todat', '-ascii')
		end
    end
	
    if options.save
        save(pwd+"\saved_data\"+dataset+".mat",'etxs','es','ens','eps','new_time_sampling','sampled_es','sampled_ens','sampled_eps' )
        dataset_id = find(dataset{:}==[datasets{:}]);
        for i=1:n_runs
            for j =1:length(methods)
                minerror(j,n_runs*(dataset_id-1)+i)=min(sampled_es{(i-1)*length(methods)+j});
            end
        end 
    end
    
	%% Displaying results
    if options.display
		set(0, 'DefaultAxesFontSize', 18);
		set(0, 'DefaultLineLineWidth', 2);

		figure; 
		for kk = 1:length(methods)
			semilogy(new_time_sampling,Means{1,kk},symbols{kk}); hold on
		end
		ylim([-inf inf])
		xlim([-inf inf])
		legend(methods_computed{:})
		xlabel('time (s)'); 
		title('mean of (||X-WH|| + lambda*logdet - min)/||X|| over time')
		grid on
		set(gca,'FontUnits','points','FontSize',18,'FontName','Times')

		figure; 
		for kk = 1:length(methods)
			semilogy(new_time_sampling,Means{2,kk},symbols{kk}); hold on
		end
		ylim([-inf inf])
		xlim([-inf inf])
		legend(methods_computed{:})
		xlabel('time (s)'); 
		title('mean of (||X-WH|| - min)/||X|| over time')
		grid on
		set(gca,'FontUnits','points','FontSize',18,'FontName','Times')

		figure; 
		for kk = 1:length(methods)
			semilogy(new_time_sampling,Means{3,kk},symbols{kk}); hold on
		end
		ylim([-inf inf])
		xlim([-inf inf])
		legend(methods_computed{:})
		xlabel('time (s)'); 
		title('mean of (lambda*logdet - min)/||X|| over time')
		grid on
		set(gca,'FontUnits','points','FontSize',18,'FontName','Times')
    end   
end

%% This is for the table with ranks
if options.save
    numalgo=length(methods);
    numprob = n_runs*length(datasets);
    for algo = 1 :  numalgo
        meanalgo(algo) = mean(minerror(algo,:)); 
        stdalgo(algo) = std(minerror(algo,:)); 
    end
    numbest = zeros(numalgo,numalgo); 

    for i = 1 : numprob
        [a,b] = sort( minerror(:,i) ); 
        c = [1:numalgo]; 
        % change order in case of equal
        for kk = 2 : numalgo
            if a(kk) <= a(kk-1)
                c(kk) = c(kk-1);
            end
        end
        for algo = 1 : numalgo
            numbest(b(algo),c(algo)) = numbest(b(algo),c(algo))+1; 
        end
    end
 
    % LaTeX Tables
    fileID = fopen(pwd+"\saved_data\"+'real_accuracy.txt','w');
    fprintf(fileID,'\\begin{center}  \n \\begin{table}[h!] \n \\begin{center} \n');
    fprintf(fileID,'Average error, standard deviation and ranking among the different runs. \n'); 
    fprintf(fileID,' \\begin{tabular}{|c|c|c|c|} \n \\hline');
    fprintf(fileID,' Algorithm &  mean $\\pm$ std & ranking  \\\\ \n \\hline \n');
    for algo = 1 : numalgo
        fprintf(fileID,replace(methods_computed{algo},"\_","-")+' & '); 
        fprintf(fileID,' $%1.3d \\pm %1.3d$ & ', meanalgo(algo) , stdalgo(algo) );  
        temp = repmat("%2.0f,",[1,numalgo]);
        temp = [temp{:}];
        temp = temp(1:end-1);
        fprintf(fileID,' ( '+string(temp)+' ) ', numbest(algo,:) ); 
        fprintf(fileID,'  \\\\ \n');
    end
    fprintf(fileID,'\\hline \\end{tabular} \n \\end{center} \n \\end{table} \n \\end{center} \n \n ');
    fclose(fileID);
end