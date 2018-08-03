# PcrBatching uses a nearest neighbor chain algorithm to 
# batch pcr operations with their nearest neighbors by extension time until
# CYCLER_COUNT groups remain. 
# These groups are marked by the maximum extension they contain
#
# then, for each thermocycler group, the operations within that group
# are clustered again with a nearest neighbor chain algorithm by temperature, 
# that refuses to make temp groups of larger than COLUMN_COUNT
# until ROW_COUNT temp groups remain.
# each row-grouping is marked with its average temperature.
# range of temp gradient is determined by the highest and lowest
# temp groups. Groups are placed into the rows that have nearest
# temperature to their average, rounding up always.

module PcrBatching
	include Utility

	COLUMN_COUNT = 12
	ROW_COUNT = 8
	CYCLER_COUNT = 4
	# THIS EXTENSION TIME IS TOO LOW?
	SEC_PER_KB_EXTENSION = 30 # 30 second, extension timer per KB for KAPA

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

			# build complete graph with edges between clusters as the absolute difference between their extension times 
			initial_graph = build_dissimilarity_matrix(pcr_operations) { |a, b| abs(a.extension_time - b.extension_time) }

			# remove all edges except those needed for mst, and then represent this graph as 
			# a min heap of edges, with extension time difference as the priority value
			# and adding the operations to the list represented as singleton clusters
			@adjacency_list = build_mst_adjacency_list(initial_graph, pcr_operations) { |pcr_op| ExtensionCluster.singleton_cluster(pcr_op) }
		end

		def combine_nearest_neighbors
		end

	end

	# A cluster of PCR operations based on the
	# nearness of their extension times
	class ExtensionCluster
		def initialize(size, min, max, mean, members)
			@group_size = size
			@min_extension = min
			@max_extension = max
			@mean_extension = mean
			@members = members
		end

		def self.singleton_cluster(pcr_operation)
			ext = pcr_operation.extension_time
			ExtensionCluster.new(1, ext, ext, ext, [pcr_operation])
		end

		def self.combine(a, b)
			combined_size = a.group_size + b.group_size
			combined_min = min(a.min_extension, b.min_extension)
			combined_max = max(a.max_extension, b.max_extension)
			combined_mean = combine_means(a.group_size, b.group_size, a.mean_extension, b.mean_extension)
			combined_members = a.members + b.members # replace with concat for in place array joining and a huge speed boost
			ExtensionCluster.new(combined_size, combined_min, combined_max, combined_mean, combined_members)
		end
	end
end

module Utility
	def max (a,b)
  		a>b ? a : b
	end

	def min (a,b)
  		a<b ? a : b
	end

	def combine_means(n1, n2, e1, e2)
		(n1*e1 + n2*e2).fdiv(n1 + n2)
	end

	# Build an adjacency matrix
	# representing the numerical difference between any 2 nodes
	# O(n^2) time complexity
	def build_dissimilarity_matrix(nodelist)
		matrix = Array.new(@clusters.size) { |i| Array.new(@clusters.size) }
		nodelist.each_with_index do |a, i|
			nodelist.each_with_index do |b, j|
				matrix[i][j] = yeild(a, b)
			end
		end
		matrix
	end

	def build_mst_adjacency_list(graph, nodelist)
		parent = prim(graph)

		adjacency_list = FastContainers::PriorityQueue.new(:min)
		(nodelist.size-1).times do |i| i + 1
			adjacency_list.push([yield(nodelist[i]), yield(nodelist[j])], graph[i][parent[i]])
		end
		adjacency_list
	end

	# naive prims algorithm to find minimum spanning tree
	# O(n^2)
	# @param graph [Array<Array<Float>>]  matrix holding the distance of every pair of nodes in the graph, 
	#                           with a distance of -1 for self-pairings
	# @return [Array<Integer>]  disjoint-set forest where each item id can be used to traverse back up the path travelled 
	def prim(graph)
		n = graph.size
		parent = []
		key = []
		visited = []

		n.times do |i|
			key[n] = Float.MAX
			visited[n] = false
		end

		key[0] = 0
		parent[0] = -1
		
		(n-1).times do
			i = min_key(key, visited)
			visited[i] = true
			n.times do |j|
				if graph[i][j] >= 0 && visited[j] == false && graph[i][j] < key[j]
					parent[j] = i
					key[j] = graph[i][j]
				end
			end
		end
		parent
	end

	# finds the index of the minimum value in an array
	def min_key key, visited
		min = Float.MAX
		mindex = -1
		key.each_with_index do |el, idx|
			if el < min && visited[idx] == false
				min = el
				mindex = idx
			end
		end
		mindex
	end

end