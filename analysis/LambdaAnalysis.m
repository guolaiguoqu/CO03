classdef LambdaAnalysis < Analysis
    properties(SetAccess=?Vertex)
        mass;  % Mass histogram
    end
    properties(SetAccess=?Helix)
        helix;
    end
    properties(Constant)
        m_pi=0.13957;  % Pion mass in Gev/c^2
        m_pro=0.9383;  % proton mass in Gev/c^2
    end
    methods
        %constructor?
        function start(obj)
            %obj.mass=Histogram(100,1.0,1.2);
            % Create mass histogram: 100 bins from 1.0 to 1.2 GeV/c^2
            obj.mass=Histogram(100,1.05,1.3);        
        end

        function event(obj,ev)
            %Author: Zichi Zhang
            %Date: 21st Feb 2026           
            % ev is a single event, obj is a LambdaAnalysis object
            %this function fill the mass and ct histograms            
            if isempty(ev) % avoid bugs caused by inconsistence with superclass
                return;
            else
                wpvxy=ev.vertex; % location of primary vertex
                tracks=ev.tracks; % get track array
                ntr=numel(tracks); % number of tracks

                % loop over unordered pairs (i<j)
                for i=1:ntr-1
                    trk1=tracks(i);

                    % Precompute per-i (cheap) quantities once:
                    q1=sign(trk1.curvature);
                    %sign of charge for the track
                    h1=Helix(trk1);
                    pt1=h1.pt();
                    % transverse momentum of the helix
                    d01=abs(h1.ipt(wpvxy));
                    % true impact parameter of the track (w.r.t. primary vertex)

                    if d01<=0.2
                        continue; % require per-track impact parameter > 0.2 cm
                    end

                    for j=i+1:ntr % avoid overcounting
                        %similiar structure with i:
                        trk2=tracks(j);
                        q2=sign(trk2.curvature);
                        if q1==q2
                            continue; % require opposite charge
                        end

                        h2=Helix(trk2);
                        d02=abs(h2.ipt(wpvxy));
                        if d02<=0.2
                            continue;
                        end

                        pt2=h2.pt();
                        if max(pt1,pt2)<=2.0
                            continue; % require the higher-pT track to exceed 2 GeV/c
                        end

                        vtx=Vertex(h1,h2); % primary vertex of helices
                        if isempty(vtx.loc)
                            continue; % no intersection
                        end

                        % requirement of lambda
                        Ld=vtx.Lxyd0pv(wpvxy);
                        Lxy=Ld(1); % transverse flight distance
                        d0Lambda=Ld(2); % impact parameter of candidate

                        if Lxy<=2.0
                            continue; 
                        end
                        if abs(d0Lambda)>=0.5
                            continue; 
                        end
                        %similiar with KOS

                        if pt1>=pt2
                            m_lambda=vtx.mass(obj.m_pro,obj.m_pi);
                        else
                            m_lambda=vtx.mass(obj.m_pi,obj.m_pro);
                        end
                        % assign masses: proton -> higher-pT track
                        obj.mass.fill(m_lambda); % fill analysis histogram

                    end % j
                end % i
            end %not empty
        end

        function stop(obj)
            % nothing to do here
        end
    end
end



