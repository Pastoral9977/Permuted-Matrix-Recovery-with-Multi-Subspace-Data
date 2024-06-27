function orth_v2(ns, rs)
    close
    
    if nargin < 2
        ns = 100:100:2000;
    end
    
    if nargin < 3
        rs = 4:4:20;
    end
    
    figure; 
    display_progress()
    for rid = 1:length(rs)
        r = rs(rid);
        norms = zeros(1,length(ns));
        normsorth = zeros(1,length(ns));
        
        subplot(1,2,1)
        hold on 
        
        for nid = 1:length(ns)
            n = ns(nid);
            sigma = 1/sqrt(n);
            A = randn(n, r)*sigma;
            [U,~,~] = svd(A);
            U = U(:,1:r);
%             norms(nid) = norm(A-U)/sqrt(r); % f-norm
            norms(nid) = norm(A-U,2)/norm(U,2);
            update_progress(nid+length(ns)*(rid-1)*2, 2*length(rs)*length(ns));
        end
        plot(ns, norms)
        hold off % 关闭绘图保持状态
        grid on
        xlim([min(ns), max(ns)]);
        xlabel('n')
        ylabel('f-norm')
        title('norm ratio of A-U over U')
        legend(cellstr(num2str(rs', 'r=%-d')))
        
        % 绘制第二个子图
        subplot(1,2,2)
        hold on % 开启绘图保持状态，使得每次绘图不会清除之前的图形
        for nid = 1:length(ns)
            n = ns(nid);
            sigma = 1/sqrt(n);
            A = randn(n, r)*sigma;
%             normsorth(nid) = norm(A'*A-eye(r), 'fro')/sqrt(r); 
            normsorth(nid) = norm(A'*A-eye(r), 2);
            update_progress(nid+length(ns)*(rid*2-1), 2*length(rs)*length(ns));
        end
        plot(ns, normsorth)
        hold off % 关闭绘图保持状态
        grid on
        xlim([min(ns), max(ns)]);
        xlabel('n')
        ylabel('f-norm ratio')
        title('norm ratio of AT*A-I over I')
        legend(cellstr(num2str(rs', 'r=%-d')))
        
    end
end
