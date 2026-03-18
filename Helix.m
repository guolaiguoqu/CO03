 classdef Helix
    % used for geometrical calculations on CdfTrack
    properties(SetAccess=?CdfTrack)
        trk; % a CdfTrack
    end
    properties(Constant)
        kpc = 0.002116;
    end
    methods
    % constructor
        function objarray = Helix(t) 
        %t is an array of CdfTrack objects
        % output an array of helix object with the track property of ith 
        % component =t(i)
            if nargin > 0  %number of argument>0
                objarray(numel(t)) = Helix; 
                % create an array of Helix objects with the correct size
                for i = 1:numel(t)
                    objarray(i).trk = t(i);
                end
            end
        end

        function r=radius(obj)          
            %find radius of circle of track of helix
            r=0.5/obj.trk.curvature;
        end

        function wcxy=center(obj)
        %   find the center of circle at POCA in xy plane
        % wcxy is the 2*1 matrix [wcx,wcy]
        % to be improved
            r=obj.radius();
            d0=obj.trk.d0;
            phi0=obj.trk.phi0;
            wc=(r+d0)*exp(1i*(phi0+pi/2)); 
            wcxy=[real(wc),imag(wc)];
            %1i to stand for imaginary number
        end
        
        function p=points(obj,t)
        % t is the angle made in xy plane 
        % between the radius connecting the point of interest
        % and the radius connecting POCA
            r=obj.radius();
            phi0=obj.trk.phi0;
            wcxy=obj.center();
            wc=wcxy(1)+1i*wcxy(2); %complex representation of center
            z0=obj.trk.z0;
            cotTheta=obj.trk.cotTheta;
            q=sign(r); %sign of charge. sign of radius. 
            w=wc-1i*r*exp(1i*(phi0+q*t));
            %size(w)=size(t)
            wx=real(w);
            wy=imag(w);
            z=z0+abs(r)*t*cotTheta;
            p=[wx; wy; z]; 
            % we need a column vector here
            % with multiple inputs (size(t)=1*n)
            %size(p)=3*n, as desired
        end
        
        function x=intersect(obj,h)
             %Author: Zichi Zhang
             %Date: 18th Feb 2026            
             %h is a helix for which we want to find the intersection with
             %obj
             %x would be a 2-column matrix giving information about 
             % two possible intersections
             %first three row contains x,y,z coordinates
             %fourth row contains difference in z at intersection
             %the fifth and sixth row contains t at the intersection
             %viewed from obj and h respectively
             if ~isa(h, 'Helix')
             %isa(object,classname) returns 1 if object is in this class
             %~ means "not"
             %if h is not a Helix object
                 disp('Argument of intersect() is not a Helix');
                 x = [ ]; % no intersection
                 return;
             end

             %whether there is intersection in xy plane:
             c1xy = obj.center();
             c1=c1xy(1)+1i*c1xy(2);
             c2xy = h.center();
             c2=c2xy(1)+1i*c2xy(2);
             d = c2xy - c1xy; %displacement vector between centers
             d2 = d(1)*d(1) + d(2)*d(2); % square of distance between centers
             r1 = obj.radius();
             r2 = h.radius();
             q1=sign(r1);
             q2=sign(r2);
             if (d2 > (abs(r1)+abs(r2))^2)
                 %equavalent to abs(d)>abs(r1)+abs(r2)
                 %but we cannot take abs for a vector
                 x = [ ]; % no intersection
                 return;
             end

             cosAlpha=(d2+r1^2-r2^2)/(2*abs(r1)*sqrt(d2));
             %cos rule
             Alpha=acos(cosAlpha);
             %angle made by the line connecting centers and 
             %the line from c1 to intersection
             Delta = atan2(c2xy(2)-c1xy(2), c2xy(1)-c1xy(1));
             %gives correct the quadrant 

             wplus=c1+abs(r1)*exp(1i*(Alpha+Delta));
             wminus=c1+abs(r1)*exp(1i*(-Alpha+Delta));
             %complex representations of intersections in xy plane

             phi0obj=obj.trk.phi0; 
             phi0h=h.trk.phi0;
             %the phi0 for each of the helix

             %use w=wc-irexp(i(phi0+qt)):
             t1plus=(angle((c1-wplus)/1i/r1)-phi0obj)*q1;
             p1plus=points(obj,t1plus);
             t1minus=(angle((c1-wminus)/1i/r1)-phi0obj)*q1;
             p1minus=points(obj,t1minus);
             t2plus=(angle((c2-wplus)/1i/r2)-phi0h)*q2;
             p2plus=points(h,t2plus);
             t2minus=(angle((c2-wminus)/1i/r2)-phi0h)*q2;
             p2minus=points(h,t2minus);
             %calculate coordinates of both points on both helices

             pplusmean=(p1plus+p2plus)/2;
             pminusmean=(p1minus+p2minus)/2;
             %averaged intersection
             dzplus=abs(p2plus(3)-p1plus(3));
             dzminus=abs(p2minus(3)-p1minus(3));
             %difference in z

             x=[[pplusmean;dzplus;t1plus;t2plus],[pminusmean;dzminus;t1minus;t2minus]];
             return;
             %t is the actual angle minus angle of POCA (with x axis)
        end

        function x = pt(obj) %momentum magnitude (defined positive)
            %unit is Gev/c
            x = obj.kpc / abs(obj.trk.curvature); %why?
            %by equating centralpetal force with magnetic force
            %(assuming uniform magnetic field in z diretion)
            %one concludes sizept is proportional to radius
        end

        function m=pvector(obj,t) %momentum vector
            %t is the angle made with POCA
            sizept=obj.pt();
            h=obj.trk.curvature;
            q=sign(h);
            cotTheta=obj.trk.cotTheta;
            phi0=obj.trk.phi0;
            %redefine variables for convenience
            px=sizept*cos(phi0+q*t);
            %phi0+q*t is the angle made by momentum vector with x axis
            py=sizept*sin(phi0+q*t);
            pz=sizept*cotTheta;
            m=[px;py;pz];
        end
        %where should this function be located?
        %should it be in vertex or helix class?

        function dpv=ipt(obj,wpvxy) %true impact parameter
            %wpvxy is the transverse location of primary vertex(true interaction)
            %primary vertex is a property of CdfEvent
            r=obj.radius();
            wcxy=obj.center();
            wc=wcxy(1)+1i*wcxy(2); 
            wpv=wpvxy(1)+1i*wpvxy(2);
            %wpv is the complex representation
            dpv=abs(wc-wpv)-abs(r);
        end
        
    end
end % classdef

           % r=obj.radius;
         %   d0=obj.trk.d0;
