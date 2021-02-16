saved_mat = ["20news" "reviews" "sports" "urban" "indian" "terrain" "sandiego" "paviau"];
root = "11-02-2021";
td = zeros(1,length(saved_mat));
i = 0;
for mat = saved_mat
    i = i+1;
    load(root+"\"+mat+".mat")
    n_runs = length(sampled_ens)/2;
    temp1 = [sampled_ens{1:2:2*n_runs}];
    temp1 = mean(temp1,2);
    temp2 = [sampled_ens{2:2:2*n_runs}];
    temp2 = mean(temp2,2);
    min2 = min(temp2);
    if min2 < min(temp1)
        td(i) = NaN;
    else
        gt0 = temp1-min2;
        gt0 = gt0>0;
        mask = [-1; 1];
        zero_crossing = conv2(gt0,mask,'same');
        idx = find(zero_crossing,1,'last');
        td(i) = new_time_sampling(end) - new_time_sampling(idx);
    end
end

% LaTeX Tables
fileID = fopen(pwd+"\saved_data\"+'timediff.txt','w');
fprintf(fileID,'\\begin{center}  \n \\begin{table}[h!] \n \\begin{center} \n');
fprintf(fileID,' \\begin{tabular}{|c|c|} \n \\hline');
fprintf(fileID,' Dataset & Time advance \\\\ \n \\hline \n');
i = 0;
for dataset = saved_mat
    i = i+1;
    fprintf(fileID, dataset + " & ");
    fprintf(fileID,' $%3.0f$ ', td(i) );  
    fprintf(fileID,'  \\\\ \n');
end
fprintf(fileID,'\\hline \\end{tabular} \n \\end{center} \n \\end{table} \n \\end{center} \n \n ');
fclose(fileID);
