classdef K0SAnalysis < Analysis
    properties(SetAccess=?Vertex)
        mass;  % Mass histogram

        ct; %lifetime of K0S signal
        ctside; %noise
        lxy;
        lxyside;
    end
    properties(SetAccess=?Helix)
        helix;
    end
    properties(Constant)
        m_pi = 0.13957;  % (charged) Pion mass in Gev/c^2
        m_K0Sworld=0.497614;
    end
    properties(Access=public) %allow changing externally
        minlxy;   % minimum Lxy 
        maxd0;    % maximum |d0| allowed 
    end

    methods
        %constructor?
        function start(obj)
            % Create mass histogram: 100 bins from 0.4 to 0.6 GeV/c^2:
            obj.mass = Histogram(100,0.4,0.6);
            obj.ct = Histogram(100,0,30); % parameter to be adjusted
            obj.ctside = Histogram(100,0,30);
            obj.lxy = Histogram(100,0,150); % parameter to be adjusted
            obj.lxyside = Histogram(100,0,150);          
            if isempty(obj.minlxy) %allow external change
                obj.minlxy = 2.0;   % default require Lxy > 2.0 cm
            end
            if isempty(obj.maxd0)
                obj.maxd0  = 0.5;   % default require |d0| < 0.5 cm
            end

        end


        function event(obj,ev)
            %Author: Zichi Zhang
            %Date: 20th Feb 2026           
            % ev is a single event, obj is a KOSAnalysis object
            %this function fill the mass and ct histograms
            if isempty(ev)
                return;
            end
            wpvxy = ev.vertex;   % location of primary vertex
            tracks = ev.tracks;  % get track array
            ntr = numel(tracks); % number of tracks

            % loop over unordered pairs (i<j):
            for i = 1:ntr-1
                %i-related conditions:
                trk1 = tracks(i);
                q1 = sign(trk1.curvature);
                h1 = Helix(trk1);
                pt1 = h1.pt();
                %magnitude of horizontal momentum
                if pt1<= 1.0
                    continue;
                end
                d01=h1.ipt(wpvxy);
                if abs(d01) <= 0.3
                    continue;
                end
                % d01 here is the true impact parameter
                % it is not just the impact parameter of track
                % loop over j for this i:
                for j = i+1:ntr
                    trk2 = tracks(j);
                    q2 = sign(trk2.curvature);
                    if q1 == q2
                        continue;
                    end
                    %  j-related cheap check first: require opposite charge 

                    h2 = Helix(trk2);
                    pt2 = h2.pt();

                    % magnitude of transverse momentum of pions (j-related)
                    if pt2<= 1.0
                        continue;
                    end

                    d02=h2.ipt(wpvxy);
                    if abs(d02) <= 0.3
                        continue;
                    end 
                    %requirement on impact parameter of track

                    vtx = Vertex(h1, h2);   % primary vertex of helices
                    if isempty(vtx.loc)
                        continue;   % no intersection
                    end

                    Ld = vtx.Lxyd0pv(wpvxy);
                    if numel(Ld) < 2
                        continue; %empty Ld for some reasons
                    end
                    Lxy = Ld(1);      % transverse flight distance
                    d0_K0S = Ld(2);   % impact parameter

                    % constraint on K0S track
                    if Lxy<=obj.minlxy
                        continue;
                    end
                    if abs(d0_K0S)>=obj.maxd0
                        continue;
                    end
                    m_K0S = vtx.mass(obj.m_pi, obj.m_pi);
                    obj.mass.fill(m_K0S);
                    % fill mass analysis histogram
     

                    if m_K0S<0.46||m_K0S>0.54
                        continue;
                    end
 
                    t1=vtx.t(1); %the angle arguments
                    t2=vtx.t(2);
                    pv1=h1.pvector(t1);
                    pv2=h2.pvector(t2);
                    pvt=pv1+pv2;
                    pT=norm(pvt(1:2)); % unit= Gev/c
                    %calculate magnitude of horizontal momentum of K0S
                    %tao_K0S=Lxy/pT*m_K0S;
                    tao_K0S=Lxy/pT*obj.m_K0Sworld;  % unit= cm/c

                    if m_K0S>=0.48&&m_K0S<=0.52
                        obj.ct.fill(tao_K0S);
                        obj.lxy.fill(Lxy);
                    else
                        obj.ctside.fill(tao_K0S);
                        obj.lxyside.fill(Lxy);
                    end
                    %lifetime histogram
                end % j
            end % i
        end

        function stop(obj)
            % nothing to do here
        end
    end
end




