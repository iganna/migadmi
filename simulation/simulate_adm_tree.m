% initial tree of n populations

path_simdata = 'sim_data/in2/';
mkdir(path_simdata)

n_leaves = 4;
names_sources = arrayfun(@(x) ['source' num2str(x)], 1:n_leaves, 'UniformOutput', false);
n_admixtures = 2;
names_admixtures = arrayfun(@(x) ['admix' num2str(x)], 1:n_admixtures, 'UniformOutput', false);
names_leaves = [names_sources, names_admixtures]
n_max = 1000;
x_thresh = log((1-1/n_max)/(1/n_max));
s = 0.1;


n_sim = 50;
f_possible = 0.01:0.01:0.99

% sf = scaling_factors
sf = [2, 1, 0.5, 0.3, 0.2];

corr_d = [];
for i_sim = 1:n_sim
    
        
    file_admix = [path_simdata 'admix', num2str(i_sim, '%02i') '.txt']
    file_len_route = [path_simdata 'len_routes', num2str(i_sim, '%02i') '.txt']
    file_alphas = [path_simdata 'alphas', num2str(i_sim, '%02i') '.txt']
    file_idxs = [path_simdata 'idxs', num2str(i_sim, '%02i') '.txt']
%     file_freqs = [path_simdata 'freqs',  num2str(i_sim, '%02i') '_' num2str(round(f_init*100), '%02i') '.txt']
    file_freqs = [path_simdata 'freqs',  num2str(i_sim, '%02i') '.txt'];
    file_d = [path_simdata 'd',  num2str(i_sim, '%02i') '.txt'];
    file_fig = [path_simdata 'fig',  num2str(i_sim, '%02i') '.pdf'];
    
    t = rand(6,1);
    len_routes = [[0;0;0;0],[t(1);t(1);t(2);t(2)] * sf(1), [t(3);t(4);t(5);t(6)]*sf(2)];
    len_routes_init = size(len_routes, 2);
    admixtures = [];
    alphas = [];
    idxs_sources = [];
    % Admixtures
    fclose(fopen(file_admix, 'w'));
    
    dmx = []
    for i = 1:4
        for j = 1:4
            idx = len_routes(i,:) ~= len_routes(j,:);
            dmx(i,j) = sum(sum(len_routes([i,j], idx)))
        end
    end

    for i_adm = 1:n_admixtures
        n_route = size(len_routes, 2);
        n_pop = size(len_routes, 1);
        
%         n_sources = randi(n_pop - 3) + 1;
        n_sources = randi(2) + 1;
        idx_sources = randperm(n_pop); idx_sources = idx_sources(1:n_sources); idx_sources = sort(idx_sources);
        adm_alpha = rand(n_sources, 1); adm_alpha = adm_alpha/sum(adm_alpha);
        
        
        %         route_new = max(max(cumsum(len_routes,2))) + 0.1;
        route_new = sum(max(cumsum(len_routes(idx_sources,:),2)')' .* adm_alpha);
        route_new = [route_new, repmat(0, 1, n_route-1)]
        len_routes = [len_routes; route_new];
        last_segment = len_routes(:, size(len_routes, 2));
        last_segment(last_segment == 0) = len_routes(last_segment == 0, 1)
        
        len_routes = [len_routes, last_segment * sf(2+i_adm)];
        
        dmx_new = mean(dmx(idx_sources, :) .* repmat(adm_alpha, 1, size(dmx, 2)));
        dmx = [dmx; dmx_new]; dmx = [dmx, [dmx_new'; 0]];
        dmx = dmx + repmat(len_routes(:,size(len_routes, 2)), 1, size(dmx, 1));
        dmx = dmx + repmat(len_routes(:,size(len_routes, 2)), 1, size(dmx, 1))';
        dmx = dmx - diag(diag(dmx));
        
       
        % Remember
        alphas = [alphas; {adm_alpha}];
        idxs_sources = [idxs_sources; {idx_sources}];
        adm_sources = names_leaves(idx_sources);
        admixtures = [admixtures; {adm_sources}];

        % Save
        fileID = fopen(file_admix,'a');
        fprintf(fileID,'%s: %s\n', names_admixtures{i_adm},...
            strjoin(adm_sources,', '));
        fclose(fileID);

    end    
    
    writematrix('', file_alphas, 'Delimiter', '\t')
    
    for itmp = 1:length(alphas)
        dlmwrite(file_alphas, alphas{itmp},'-append', 'Delimiter', '\t')    
    end
    
    
    
    writecell(idxs_sources, file_idxs)
    
    

    
    % Points to visualize only, not for frequencies
    points = [[0;0;0;0],[1.5;1.5;-1.5;-1.5],[2;1;-1;-2]];
    for i_adm = 1:n_admixtures
        idx_sources = idxs_sources{i_adm};
        adm_alpha = alphas{i_adm};
        p_sources = points(idx_sources, size(points, 2))
        p_new = sum(adm_alpha .* p_sources)
        p_new = [repmat(p_new, 1, size(points, 2))]
        points = [points; p_new];
%         points = [points, points(:,size(points, 2)) + (rand(size(points, 1), 1) - 0.5)*0.1];
        points = [points, points(:,size(points, 2))];
    end
    
    
    
    if(i_sim == 1)
        % Figure
        len_routes_cum = cumsum(len_routes,2);
        f1 = figure; hold on;
        for i = 1:size(points, 1)
            h = plot(len_routes_cum(i,:), points(i,:), '>-', 'LineWidth', 1.5);
            h.MarkerFaceColor = h.Color;
        end
        for i_adm = 1:n_admixtures
            idx_sources = idxs_sources{i_adm};
            i_pos = len_routes_init + i_adm - 1;
            for j_adm = idx_sources
                x = [len_routes_cum(j_adm,i_pos), len_routes_cum(i_adm + n_leaves,i_pos)];
                y = [points(j_adm,i_pos), points(i_adm + n_leaves,i_pos)];
                h = quiver(x(1),y(1),x(2) - x(1),y(2) - y(1),0.97,'Color', 'k', ...
                    'MaxHeadSize',0.3/sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2), 'LineStyle', '-')
            end
        end
        legend(names_leaves)
        saveas(gcf, file_fig);
        close(f1);
% 
    end
    
    writematrix(len_routes, file_len_route)
        
    % Freqs
    f_all = [];
    for f_init = f_possible
%     for f_init = [0.2, 0.3, 0.5, 0.7, 0.8]
            
        x_init = log((1-f_init)/f_init);
        f_branch = repmat(f_init, n_leaves, 1);
        x_branch = repmat(x_init, n_leaves, 1);

        for i = 2:len_routes_init


            [v_un, ia, ib] = unique(len_routes(1:n_leaves,i));
            f_un = f_branch(ia,i-1);
            x_un = x_branch(ia,i-1);
            v_mul = max((f_init .* (1-f_init)), (f_un .* (1-f_un)));
            x_add = normrnd(0, s * v_un ./ v_mul);

            x_new = x_branch(:,i-1) +  x_add(ib);
            x_new = max(min(x_new, x_thresh), -x_thresh);
            x_new(x_branch(:,i-1) == x_thresh) = x_thresh;

            x_branch(:,i) = x_new;
            f_branch(:,i) =  1 ./ (exp(x_branch(:,i)) + 1);

            max(abs(f_branch(:,i) - f_branch(:,i-1))); 

            if max(abs(f_branch(:,i) - f_branch(:,i-1))) > 0.5
                break
            end
        end


        for i_adm = 1:n_admixtures
            % Get admixture parameters
            idx_sources = idxs_sources{i_adm};
            adm_alpha = alphas{i_adm};
            % Generate new value
            f_sources = f_branch(idx_sources, size(f_branch, 2))
            x_sources = log((1-f_sources)./f_sources)
            x_new = sum(adm_alpha .* x_sources)
            f_new = 1./ (exp(x_new) + 1)
            f_new = [repmat(f_new, 1, size(f_branch, 2))];
            % Update
            f_branch = [f_branch; f_new];
            
            
            v_un = len_routes(1:size(f_branch,1),size(f_branch, 2)+1);
            f_un = f_branch(:,size(f_branch, 2));
            x_un = log((1-f_un)./f_un);
            v_mul = max((f_init .* (1-f_init)), (f_un .* (1-f_un)));
            x_add = normrnd(0, s * v_un ./ v_mul);

            x_new = x_un +  x_add;
            x_new = max(min(x_new, x_thresh), -x_thresh);
%             x_new(x_branch(:,i-1) == x_thresh) = x_thresh;

            
            f_branch_new =  1 ./ (exp(x_new) + 1);
            
            f_branch = [f_branch, f_branch_new];
        end
        

%         len_routes_cum = cumsum(len_routes,2);
%         f1 = figure; hold on;
%         for i = 1:size(f_branch, 1)
%             plot(len_routes_cum(i,:), f_branch(i,:), 'o-');
%         end
%         ylim([0 1]);

        % saveas(gcf, [path_freqs_tmp 'sim.pdf']);
        % close(f1);
        f_all = [f_all; f_branch(:,size(f_branch, 2))'];
    end
    
    writecell(names_leaves, file_freqs, 'Delimiter', '\t')
    dlmwrite(file_freqs,f_all,'-append', 'Delimiter', '\t')
    
%     d = array2table(squareform(pdist(f_all')))
    
    d = array2table(squareform(pdist(log((1-f_all') ./f_all') )));
%     d = array2table(dmx)

    d.Properties.RowNames = names_leaves;
    d.Properties.VariableNames = names_leaves;
    writetable(d, file_d, 'Delimiter', '\t', 'WriteRowNames', true)
    
    
    tmp = table2array(d);
    tmp = tmp(:);
    corr_d = [corr_d; corr(tmp(:), dmx(:))];
    
end

% 
% d = squareform(pdist(f_all'))
% 
% t = seqneighjoin(d(1:n_leaves, 1:n_leaves),'equivar', names_sources)
% plot(t)
% 
% t = seqneighjoin(d,'equivar', names_leaves)
% plot(t)



