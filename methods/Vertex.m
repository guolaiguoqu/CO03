classdef Vertex             
    %note this class is about the intersection of daughter particle helices
    %not to be confused with primary vertex (about where the particle is formed)
    properties(SetAccess=?Helix) %because we need to use functions on helix
        helix; %helices 
        loc; %vertex location in xy  2*1
        dz; %difference in z of the two helices (maybe not required)
        t; % angle at the vertex w.r.t both vertices  2*1
    end
    properties(Constant)
        kpc = 0.002116;
    end
    methods
        %contructor
        function obj=Vertex(h1,h2) 
            %h1 and h2 are two helix objects
            %returns the position vector with deltaz

            inter=intersect(h1,h2);
            if ~isempty(inter) %check intersection exist
                if inter(4,1)<inter(4,2) %compare dz
                    vertex=inter(:,1); %vertex
                else
                    vertex=inter(:,2);
                end            
                obj.helix=[h1,h2];
                obj.loc=vertex(1:3);
                obj.dz=vertex(4);
                obj.t=vertex(5:6);
            else
                obj.helix=[];
                obj.loc=[];
                obj.t=[];
            end
            return;
        end
        
        function mt=mass(obj,m1,m2)
        %x is the mass of a parent particle that decayed into two tracks.
        %all the masses given are rest mass
        %two mass inputs are required because the identity of decayed
        %particles is changable
        %all quantities below use natural unit
        h1=obj.helix(1);
        h2=obj.helix(2);
        t1=obj.t(1);
        t2=obj.t(2);
        p1v=h1.pvector(t1); 
        p2v=h2.pvector(t2);
        p1=norm(p1v);
        p2=norm(p2v);
        ptv=p1v+p2v; % total momentum vector
        pt=norm(ptv);
        %evaluate magnitude of individual and total momentum

        E1=sqrt(p1^2+m1^2);
        E2=sqrt(p2^2+m2^2);
        Et=E1+E2;
        %calculate energies assuming natural units

        mt=sqrt(Et^2-pt^2);
        end

        function x=Lxyd0pv(obj,wpv) 
            %find transverse flight distance Lxy 
            %and impact parameter of primary vertex d0pv
            %note this d0pv is not d0 in CdfTrack
            %output x=[Lxy,d0pv] 
            h1=obj.helix(1);
            h2=obj.helix(2);
            t1=obj.t(1);
            t2=obj.t(2);
            p1v=h1.pvector(t1); 
            p2v=h2.pvector(t2);
            ptv=p1v+p2v; % total momentum vector
            %repetition of code: possibility of improvement?
            pT=ptv(1:2); 
            %note ptv(1,2) is the single element in 1st row 2nd column
            pT_unit=pT/norm(pT);
            loc_xy=obj.loc(1:2); %location in xy plane
            deltaloc=loc_xy-wpv;
            Lxy=dot(deltaloc,pT_unit);
            d0pv=deltaloc(1)*pT_unit(2)-pT_unit(1)*deltaloc(2); 
            %magnitude of cross product
            x=[Lxy,d0pv];
        end

    end
end

