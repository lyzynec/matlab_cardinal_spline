classdef CardinalSpline
    %TRAJECTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        nodes
        tangents
        distances
    end
    
    methods
        function obj = CardinalSpline(nodes, scale, dist_step_acc)
            %TRAJECTORY Construct an instance of this class
            %   Detailed explanation goes here

            assert(size(nodes, 1) > 2, "Argument nodes must have size of n*2 where n >= 2.");
            assert(size(nodes, 2) == 2, "Argument nodes must have size of n*2.");
            assert(isscalar(scale), "Argument scale must be scalar.");
            

            obj.nodes = nodes;
            obj.tangents = obj.calc_tangents(nodes)*scale;
            obj = obj.calc_distances(dist_step_acc);
        end


        function positions = get_positions(obj,t)
            assert(isvector(t));
            positions = zeros(length(t), size(obj.nodes,2));
            t = t';

            indexes = floor(t)+1;
            lt = mod(t,1);
            positions(lt==0, :) = obj.nodes(indexes(lt==0),:);
            M = [1,0,0,0;0,1,0,0;-3,-2,3,-1;2,1,-2,1];
            lt_matrix = reshape([ones(sum(lt~=0), 1), lt(lt ~= 0), lt(lt ~= 0).^2, lt(lt ~= 0).^3]', 1, 4, []);
            p = [
                reshape(obj.nodes(indexes(lt ~= 0),:)', 1,size(obj.nodes,2),[]);
                reshape(obj.tangents(indexes(lt ~= 0),:)', 1,size(obj.nodes,2),[]);
                reshape(obj.nodes(indexes(lt ~= 0)+1,:)', 1,size(obj.nodes,2),[]);
                reshape(obj.tangents(indexes(lt ~= 0)+1,:)', 1,size(obj.nodes,2)',[]);
            ];
            positions(lt~=0, :) = reshape(permute(pagemtimes(pagemtimes(lt_matrix, M), p), [3,2,1]), sum(lt~=0), 2);
        end

        function distances = get_distances(obj)
            distances = obj.distances;
        end

        function distance = get_distance_sum(obj)
            distance = sum(obj.distances);
        end

        function tangents = calc_tangents(~, nodes)
            len_n = size(nodes,1);
            tangents = zeros(size(nodes));
            tangents(1, :) = 2*(nodes(2, :)-nodes(1,:));
            tangents(len_n, :) = 2*(nodes(len_n,:)-nodes(len_n-1,:));
            i = 2:len_n-1;
            tangents(i, :) = nodes(i+1,:) - nodes(i-1,:);
        end

        function obj = calc_distances(obj, dist_step_acc)
            len_nodes = size(obj.nodes, 1);
            len_distances = len_nodes - 1;
            

            obj.distances = zeros(len_distances, 1);
            for distance_index = 1:len_distances
                prew_distance = 0;
                num_of_points = 2;
                while true
                    l = linspace(distance_index-1,distance_index,num_of_points);
                    positions = obj.get_positions(l);
                    
                    distance = sum_position_distances(positions);

                    if abs(distance - prew_distance) < dist_step_acc
                        obj.distances(distance_index) = distance;
                        break
                    end
                    prew_distance = distance;
                    num_of_points = num_of_points*2;

                end
            

            end

        end
    end
end

function sum_distances = sum_position_distances(positions)
    v = diff(positions, 1, 1);
    sum_distances = sum((sum(v.^2,2).^(1/2)));
end
