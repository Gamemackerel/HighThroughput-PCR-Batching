needs 'GradientPcrHelpers'

module GradientPcrRepresentation
	include GradientPcrHelpers

	# Object representation for an individual pcr reaction.
	# PcrOperation keeps track of all the necessary factors
	# for a pcr.
	class PcrOperation

		# make a brand new PcrOperation.
		#
		# @param opts [Hash]  initialization options
		# @option opts [Float] :extension_time  extension time for this pcr operation
		# @option opts [Float] :anneal_temp  the annealing temperature for this pcr operation
		# @option opts [Integer] :extension_group  a group id for this operation, shared with 
		# 							other pcr operations who could be run together in a reaction 
		# 							with the same extension time
		# @option opts [Integer] :tanneal_group  a group id for this operation, shared with 
		# 							other pcr operations who could be run together in a reaction 
		# 							with the same annealing temperature
		def initialize opts
			@extension_time 	= opts[:extension_time]
			@anneal_temp 		= opts[:anneal_temp]
			@extension_group 	= opts[:extension_group]
			@tanneal_group 		= opts[:tanneal_group]
		end

		# get an exact copy of this pcr operation
		def clone
			PcrOperation.new({
				extension_time: 	@extension_time
				anneal_temp: 		@anneal_temp
				extension_group: 	@extension_group
				tanneal_group: 		@tanneal_group
			})
		end
	end


	# A set of clusters of pcr_operations.
	# the clusters will be made by the proximity of extension
	# time, so that multiple pcr_operations can be optimally
	# put into the same pcr reaction if they have similar enough
	# extension time 
	class ExtensionClusterGraph

		# Use a list of pcr_operations to create a set of singleton clusters
		# ready for combining into larger clusters based on similarity of extension time
		#
		# @param pcr_operations [Array<PcrOperation>]  list of pcr_operations to be clustered 
		def initialize(pcr_operations)
			@size = pcr_operations.size

			# build complete graph (as adjacency matrix) with edges between 
			# clusters as the absolute difference between those clusters' extension times 
			initial_graph = build_dissimilarity_matrix(pcr_operations) { |a, b| distance_func(a,b) }

			# remove all edges except those needed for mst, and then represent this graph as 
			# a min heap of edges, with extension time difference as the priority value
			# and adding the operations to the list represented as singleton clusters
			@adjacency_list = build_mst_adjacency_list(initial_graph, pcr_operations) { |pcr_op| ExtensionCluster.singleton_cluster(pcr_op) }
		end

		def combine_nearest_clusters
			distance = @adjacency_list.peek().top_key()
			cluster_a, cluster_b = @adjacency_list.top()
			@adjacency_list.pop
			cluster_ab = ExtensionCluster.combine(cluster_a, cluster_b, distance)

			# go through adjacency list updating pairs and distances to reflect this new merge
			# lots of edge cases here ex:
			# c - a, c - b
			# after merge, c - ab, c - ab
			duplicate_checker = Set.new # two arrays with the same contents will have the same hash code
			@adjacency_list.each do |pair, priority|
				replace_index = pair.index(cluster_a) || pair.index(cluster_b)
				other_index = replace_index == 1 ? 0 : 1
				if !replace_index.nil?
					new_pair = Array.new(pair) #Priority queue probably uses hash code of the object, which is not retained for arrays on content change, so we cannot 'update this pair in the queue using its reference' 
					new_pair[replace_index] = cluster_ab 
					if duplicate_checker.contains?(new_pair) # edge case: merge will cause ab - ab pair 
						remove_heap_element(@adjacency_list, pair) #TODO DEFINE THIS with decreasekey and pop (constant)
					else
						new_priority = distance_func(pair[replace_index], pair[other_index])
						reinsert(@adjacency_list, new_pair, pair, new_priority, priority) #TODO DEFINE THIS remove_heap_element, and push (logn)
						duplicate_checker << pair
					end
				else
					duplicate_checker << pair
				end

		end



		def distance_func
			if (cluster_a.size + cluster_b.size) > (thermocycler_rows * thermocycler_columns) && (TannealCluster.anneal_range(cluster_a, cluster_b) > 10)
				# prevent combination if it would produce an anneal range or batch size that a single thermocycler cannot handle
				return Float.MAX
			else
				return abs(cluster_a.mean_extension_time - cluster_b.mean_extension_time)
			end
		end

		# decides whether or not further clustering is required
		#
		# @return [Boolean]  whether clustering has finished
		def threshhold_func force_combination_distance
			if graph_representation <= thermocycler_quantity
				if @adjacency_list.peek().priority < force_combination_distance
					false
				else
					true
				end
			else
				false
			end
		end


		def cluster_set
			clusters = Set.new
			@adjacency_list.each do |cluster_tuple, priority|
				cluster_tuple.each do |cluster|
					clusters << ExtensionCluster.get_top_level_cluster(cluster)
				end
			end
			clusters
		end
	end

	# A cluster of PCR operations based on the
	# nearness of their extension times
	class ExtensionCluster
		def initialize(opts)
			@group_size 	= opts[:size]
			@min_extension 	= opts[:min_extension]
			@max_extension 	= opts[:max_extension]
			@mean_extension = opts[:mean_extension]
			@max_anneal 	= opts[:max_anneal]
			@min_anneal 	= opts[:min_anneal]
			@members 		= opts[:member_list]
			@child_clusters = opts[:child_clusters]
			@parent_cluster = opts[:parent_cluster]
			@combo_distance = opts[:combination_distance]
		end

		def self.singleton_cluster(pcr_operation)
			ext = pcr_operation.extension_time
			min_anneal = min(pcr_operation.fwd_anneal, pcr_operation.rev_anneal)
			max_anneal = max(pcr_operation.fwd_anneal, pcr_operation.rev_anneal)
			ExtensionCluster.new(
					size: 			1, 
					min_extension: 	ext, 
					max_extension: 	ext, 
					mean_extension: ext,
					max_anneal: 	max_anneal,
					min_anneal: 	min_anneal,
					member_list: 	[pcr_operation]
				)
		end

		def self.combine(a, b, combination_distance)
			combined_size = a.group_size + b.group_size
			combined_min = min(a.min_extension, b.min_extension)
			combined_max = max(a.max_extension, b.max_extension)
			combined_mean = combine_means(a.group_size, b.group_size, a.mean_extension, b.mean_extension)
			combined_members = a.members + b.members # this is a bottleneck. 
								# replace with concat for in place array joining and a huge speed boost
								# or don't store members list, and recurse to bottom of tree to get it when needed
			ab = ExtensionCluster.new(
					size: 			combined_size, 
					min_extension: 	combined_min, 
					max_extension: 	combined_max, 
					mean_extension: combined_mean,
					max_anneal: 	max(a.max_anneal, b.max_anneal),
					min_anneal: 	min(a.min_anneal, b.min_anneal),
					member_list: 	combined_members
					child_clusters: [a,b]
					combo_distance: combination_distance
				)
			a.parent_cluster = ab
			b.parent_cluster = ab

			ab
		end

		def self.get_top_level_cluster(cluster)
			if cluster.parent_cluster.nil?
				return cluster
			else
				return get_top_level_cluster(cluster.parent_cluster)
			end
		end

		# calculate members when needed
		# lazy approach, if we don't want to keep track of the member_list for each cluster 
		def self.members(cluster)
			if child_clusters.nil?
				return [cluster.pcr_operation]
			else
				return members(cluster.child_clusters[0]).concat(members(cluster.child_clusters[1]))
			end
		end

		def to_string()
			"#{@members} \n extension range: #{min_extension}-#{max_extension} \n anneal range: #{min_anneal}-#{max_anneal} \n"
		end
	end
end