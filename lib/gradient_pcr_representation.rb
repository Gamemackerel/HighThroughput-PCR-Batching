# This file defines objects used to represent individual Pcr Operations, Clusters of Pcr Representations,
# and Graphs of Pcr operation Clusters
#

require 'PCR Libs/GradientPcrHelpers'


# Object representation for an individual pcr reaction.
# PcrOperation keeps track of all the necessary factors
# for a pcr.
#
class PcrOperation

    attr_reader :extension_time, :anneal_temp, :extension_group, :tanneal_group, :unique_id

    # make a brand new PcrOperation.
    #
    # @param opts [Hash]  initialization options
    # @option opts [Float] :extension_time  extension time for this pcr operation
    # @option opts [Float] :anneal_temp  the annealing temperature for this pcr operation
    # @option opts [Integer] :extension_group  a group id for this operation, shared with 
    #                           other pcr operations who could be run together in a reaction 
    #                           with the same extension time
    # @option opts [Integer] :tanneal_group  a group id for this operation, shared with 
    #                           other pcr operations who could be run together in a reaction 
    #                           with the same annealing temperature
    # @option opts [Integer] :unique_id  identifier to track a pcr operation
    #                           through its being batched
    def initialize opts = {}
        @extension_time     = opts[:extension_time]
        @anneal_temp        = opts[:anneal_temp]
        @extension_group    = opts[:extension_group]
        @tanneal_group      = opts[:tanneal_group]
        @unique_id          = opts[:unique_id]
        
        if @extension_time.nil? || @anneal_temp.nil? || @extension_time.zero? || @anneal_temp.zero?
            raise "Every pcr operation must have a non-zero extension time and anneal temp"
        end
    end

    # get an exact copy of this pcr operation
    def clone
        PcrOperation.new({
            extension_time:     @extension_time,
            anneal_temp:        @anneal_temp,
            extension_group:    @extension_group,
            tanneal_group:      @tanneal_group,
            unique_id:          @unique_id,
        })
    end

    def to_s 
        "extension_time: #{@extension_time} + \n anneal_temp: #{@anneal_temp} + \n extension_group: #{@extension_group} + \n tanneal_group: #{@tanneal_group}"
    end
end

# Core clustering logic that is used for both ExtensionClusterGraph and TannealClusterGraph to group
# a set of objects based on the relative proximity of their fields
#
#   representation invariants of cluster graphs:
#   initial size of pcr_operations == clusters.map { members }.flatten.size
#   Set(pcr_operations) == Set(clusters.map { members }.flatten)
#   pcr operations belong to exactly 1 cluster
#
module ClusterGraphMethods
    # combines nearest clusters until threshold function is triggered
    # modifies the state of the graph, and returns the resulting set of clusters
    def perform_clustering
        while !threshhold_func
            combine_nearest_clusters_lazy
        end
        cluster_set
    end

    # Remove the shortest distance edge and combine the cluster pair into one cluster
    # if one or both of the cluster pair has already been combined, update the edge entry 
    # to reflect the highest level cluster that they now constitute
    def combine_nearest_clusters_lazy
        distance = @adjacency_list.min_priority
        pair = @adjacency_list.delete_min_return_key
        cluster_a, cluster_b = pair.to_a

        if cluster_a.child_cluster || cluster_b.child_cluster 
            # at least one of these clusters have been merged 
            # so we must update this adjacency list entry without combining them
            a_super = cluster_a.get_containing_supercluster
            b_super = cluster_b.get_containing_supercluster
            if a_super != b_super
                @adjacency_list.push(Set[a_super, b_super], distance_func(a_super, b_super))
            end
        else
            @size -= 1
            cluster_ab = cluster_a.combine_with(cluster_b)
        end

        if @adjacency_list.empty?
            @final_cluster = cluster_ab
        end
    end

    # When using lazy combination, the adjacency list graph representation may have duplicate clusters,
    # or clusters which have already been swallowed by a larger one. `cluster_set` returns only the 
    # most current, highest level graph nodes 
    def cluster_set
        clusters = Set.new
        clusters << @final_cluster if @final_cluster
        @adjacency_list.each do |cluster_tuple, priority|
            cluster_tuple.each do |cluster|
                clusters << cluster.get_containing_supercluster
            end
        end
        clusters
    end
    
    # Gets Float value which represents an impossible combination distance.
    # Float::MAX cannot be used for this value, since it is critical to prims
    def max_distance
        (Float::MAX - 1)
    end

    def checkrep
        clusters = cluster_set
        total_pcr_members = clusters.to_a.map { |c| c.members }.flatten.size
        assert(total_pcr_members == @initial_size)
    end

    def to_s 
        "size:" + @size.to_s + "\n" + "clusters: " + cluster_set.to_a.to_s
    end
end

# Methods that are useful for both Tanneal clusters or Extension time clusters
#
module ClusterMethods
    # calculate members when needed
    # lazy approach, so we dont have to keep track of the member_list for each cluster 
    def members
        if @parent_clusters.nil?
            return [@pcr_operation]
        else
            return @parent_clusters[0].members.concat(@parent_clusters[1].members)
        end
    end

    def combine_anneal_range_with(other)
        combine_range(self.max_anneal, self.min_anneal, other.max_anneal, other.min_anneal)
    end

    def get_containing_supercluster
        if @child_cluster.nil?
            return self
        else
            return @child_cluster.get_containing_supercluster
        end
    end
end

# A set of clusters of pcr_operations.
# the clusters will be made by the proximity of extension
# time, so that multiple pcr_operations can be optimally
# put into the same pcr reaction if they have similar enough
# extension time 
#
class ExtensionClusterGraph
    include GradientPcrHelpers
    include ClusterGraphMethods
    
    attr_reader :size, :initial_size, :adjacency_list

    def initialize(opts = {})
        pcr_operations              = opts[:pcr_operations]
        @thermocycler_quantity      = opts[:thermocycler_quantity]
        @thermocycler_rows          = opts[:thermocycler_rows]
        @thermocycler_columns       = opts[:thermocycler_columns]
        @thermocycler_temp_range    = opts[:thermocycler_temp_range]
        @force_combination_distance = opts[:force_combination_distance]
        @prevent_combination_distance = opts[:prevent_combination_distance]
        @size = pcr_operations.size
        @initial_size = @size # initial size recorded for checkrep


        singleton_clusters = pcr_operations.map { |pcr_op| ExtensionCluster.singleton_cluster(pcr_op) }

        # final cluster field only stores a cluster if there is only one cluster in the graph (adjacency list cannot represent this state)
        @final_cluster = singleton_clusters.first if singleton_clusters.one?

        # build complete graph (as adjacency matrix) with edges between 
        # clusters as the absolute difference between those clusters' extension times 
        initial_graph = build_dissimilarity_matrix(singleton_clusters) do |a, b| #O(n^2)
            distance_func(a,b) 
        end

        @adjacency_list = build_mst_adjacency_list(initial_graph, singleton_clusters)  #O(n^2)
    end

    def distance_func(cluster_a, cluster_b)
        if (cluster_a.size + cluster_b.size) > (@thermocycler_rows * @thermocycler_columns)  \
           || (cluster_a.combine_anneal_range_with(cluster_b) > @thermocycler_temp_range)
            # prevent combination of pairs if it would produce an anneal range or batch size that a single thermocycler cannot handle
            return max_distance
        else
            return cluster_a.combine_extension_range_with(cluster_b)
        end
    end

    # decides whether or not further clustering is required
    #
    # @return [Boolean]  whether clustering has finished
    def threshhold_func 

        if @adjacency_list.empty?  #this case is reached when weve combined into only a single cluster
            return true
        end

        next_distance = @adjacency_list.min_priority
        next_pair = @adjacency_list.min_key
        cluster_a, cluster_b = next_pair.to_a

        # End clustering if there are no more clusters to combine, 
        # or the next combination distance is greater than or equal to 
        # the specified maximum allowable combination distance. 
        #
        # Impossible combination pairs (with distance Float.MAX) will not exceed
        # distance specified in @prevent_combination_distance
        if (next_distance >= @prevent_combination_distance)
            return true
        end

        # If we have already combined enough so that all operations can be 
        # run at once in the amount of thermocyclers available, 
        # then we combine only if the next distance is less than or equal to
        #  the specified distance for mandatory combination
        if @size <= @thermocycler_quantity
            if next_distance <= @force_combination_distance
                false
            else
                true
            end
        else
            false
        end
    end
end

# A cluster of PCR operations based on the
# nearness of their extension times
#
class ExtensionCluster
    include GradientPcrHelpers
    include ClusterMethods

    attr_reader :size, :min_extension, :max_extension, :mean_extension, :max_anneal, :min_anneal, :parent_clusters, :child_cluster, :pcr_operation
    attr_writer :child_cluster

    def initialize(opts)
        @size    = opts[:size]
        @min_extension   = opts[:min_extension]
        @max_extension   = opts[:max_extension]
        @mean_extension  = opts[:mean_extension]
        @max_anneal      = opts[:max_anneal]
        @min_anneal      = opts[:min_anneal]
        @parent_clusters = opts[:parent_clusters]
        @child_cluster   = opts[:child_cluster]
        @pcr_operation   = opts[:pcr_operation]
    end

    def self.singleton_cluster(pcr_operation)
        ext = pcr_operation.extension_time
        anneal = pcr_operation.anneal_temp
        ExtensionCluster.new(
                size:           1, 
                min_extension:  ext, 
                max_extension:  ext, 
                mean_extension: ext,
                max_anneal:     anneal,
                min_anneal:     anneal,
                pcr_operation:  pcr_operation
            )
    end

    def combine_extension_range_with(other)
        combine_range(self.max_extension, self.min_extension, other.max_extension, other.min_extension)
    end

    def combine_with(other)
        combined_size = self.size + other.size
        combined_min = min(self.min_extension, other.min_extension)
        combined_max = max(self.max_extension, other.max_extension)
        combined_mean = combine_means(self.size, other.size, self.mean_extension, other.mean_extension)
        super_cluster = ExtensionCluster.new(
                size:            combined_size, 
                min_extension:   combined_min, 
                max_extension:   combined_max, 
                mean_extension:  combined_mean,
                max_anneal:      max(self.max_anneal, other.max_anneal),
                min_anneal:      min(self.min_anneal, other.min_anneal),
                parent_clusters: [self,other]
            )
        self.child_cluster = super_cluster
        other.child_cluster = super_cluster

        super_cluster
    end

    # calculate members when needed
    # lazy approach, so we dont have to keep track of the member_list for each cluster
    # TODO: make dynamic so members only recursively calculates the first time, and then stores
    # its members list for future calls
    def members
        if @parent_clusters.nil?
            return [@pcr_operation]
        else
            return @parent_clusters[0].members.concat(@parent_clusters[1].members)
        end
    end

    def get_containing_supercluster
        if @child_cluster.nil?
            return self
        else
            return @child_cluster.get_containing_supercluster
        end
    end

    def to_s
        "size: #{@size} \n extension range: #{min_extension}-#{max_extension} \n anneal range: #{min_anneal}-#{max_anneal} \n"
    end
end

# A set of clusters of pcr_operations.
# the clusters will be made by the proximity of Annealling temperature
# so that multiple pcr_operations can be optimally
# put into the same pcr thermocycler row
# if they have similar enough annealling temperature 
#
class TannealClusterGraph
    include GradientPcrHelpers
    include ClusterGraphMethods

    attr_reader :size, :initial_size, :adjacency_list

    def initialize(opts = {})
        pcr_operations              = opts[:pcr_operations]
        @thermocycler_columns       = opts[:thermocycler_columns]
        @thermocycler_rows          = opts[:thermocycler_rows]          
        @thermocycler_temp_range    = opts[:thermocycler_temp_range]
        @force_combination_distance = opts[:force_combination_distance]
        @prevent_combination_distance = opts[:prevent_combination_distance]         
        @size           = pcr_operations.size
        @initial_size   = @size
        @final_cluster  = TannealCluster.singleton_cluster(pcr_operations.first) if pcr_operations.one?

        # build complete graph (as adjacency matrix) with edges between 
        # clusters as the absolute difference between those clusters' Tanneal 
        initial_graph = build_dissimilarity_matrix(pcr_operations) do |a, b| #O(n^2)
            distance_func(TannealCluster.singleton_cluster(a),TannealCluster.singleton_cluster(b)) 
        end

        # remove all edges except those needed for mst, and then represent this graph as 
        # a min heap of edges, with Tanneal time difference as the priority value
        # and adding the operations to the list represented as singleton clusters
        singleton_clusters = pcr_operations.map { |pcr_op| TannealCluster.singleton_cluster(pcr_op) }
        @adjacency_list = build_mst_adjacency_list(initial_graph, singleton_clusters)  #O(n^2)
    end

    def distance_func(cluster_a, cluster_b)
        if (cluster_a.size + cluster_b.size) > (@thermocycler_columns)
            # prevent combination if it would produce a tanneal group too big to fit in one row
            return max_distance
        else
            return cluster_a.combine_anneal_range_with(cluster_b)
        end
    end

    # decides whether or not further clustering is required
    #
    # @return [Boolean]  whether clustering has finished
    def threshhold_func
        next_distance = @adjacency_list.min_priority
        next_pair = @adjacency_list.min_key
        cluster_a, cluster_b = next_pair.to_a

        # End clustering if there are no more clusters to combine, 
        # or the next combination distance is greater than or equal to 
        # the specified maximum allowable combination distance. 
        #
        # Impossible combination pairs (with distance Float.MAX) will not exceed
        # distance specified in @prevent_combination_distance
        if @adjacency_list.empty? || next_distance >= @prevent_combination_distance
            return true
        end

        # If we have already combined enough so that all operations can be 
        # run at once in the amount of rows available, 
        # then we combine only if the next distance is less than or equal to
        # the specified distance for mandatory combination
        if @size <= @thermocycler_rows
            if next_distance <= @force_combination_distance
                false
            else
                true
            end
        else
            false
        end
    end
end

# A cluster of PCR operations based on the
# nearness of their Tanneal times
#
class TannealCluster
    include GradientPcrHelpers
    include ClusterMethods

    attr_reader :size, :min_anneal, :max_anneal, :mean_anneal, :parent_clusters, :child_cluster, :pcr_operation
    attr_writer :child_cluster

    def initialize(opts)
        @size            = opts[:size]
        @min_anneal      = opts[:min_anneal]
        @max_anneal      = opts[:max_anneal]
        @mean_anneal     = opts[:mean_anneal]
        @parent_clusters = opts[:parent_clusters]
        @child_cluster   = opts[:child_cluster]
        @pcr_operation   = opts[:pcr_operation]
    end

    def self.singleton_cluster(pcr_operation)
        anneal = pcr_operation.anneal_temp
        TannealCluster.new(
                size:           1, 
                min_anneal:     anneal, 
                max_anneal:     anneal, 
                mean_anneal:    anneal,
                pcr_operation:  pcr_operation
            )
    end

    def combine_with(other)
        combined_size = self.size + other.size
        combined_min = min(self.min_anneal, other.min_anneal)
        combined_max = max(self.max_anneal, other.max_anneal)
        combined_mean = combine_means(self.size, other.size, self.mean_anneal, other.mean_anneal)
        super_cluster = TannealCluster.new(
                size:            combined_size, 
                min_anneal:      combined_min, 
                max_anneal:      combined_max, 
                mean_anneal:     combined_mean,
                parent_clusters: [self,other]
            )
        self.child_cluster = super_cluster
        other.child_cluster = super_cluster

        super_cluster
    end

    def to_s
        "size: #{@size} \n Tanneal range: #{@min_anneal}-#{@max_anneal}"
    end
end