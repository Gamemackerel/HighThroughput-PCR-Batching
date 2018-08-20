module GradientPcrHelpers

	require 'priority_queue'

	def assert(expr)
		raise "This is wrong" unless expr
	end

	def max (a,b)
  		a>b ? a : b
	end

	def min (a,b)
  		a<b ? a : b
	end

	def combine_means(n1, n2, e1, e2)
		(n1*e1 + n2*e2).fdiv(n1 + n2)
	end


	def combine_range(a_max, a_min, b_max, b_min)
		max(a_max - b_min, b_max - a_min)
	end

	# Build an adjacency matrix
	# representing the numerical difference between any 2 nodes
	# accepts code block for |a,b| difference function
	# O(n^2) time complexity
	def build_dissimilarity_matrix(nodelist)
		matrix = Array.new(nodelist.size) { |i| Array.new(nodelist.size) }
		nodelist.each_with_index do |a, i|
			nodelist.each_with_index do |b, j|
				matrix[i][j] = yield(a, b)
			end
		end
		matrix
	end

	# remove all edges except those needed for mst, and then represent this graph as 
	# a min heap of edges, with extension time difference as the priority value
	# and adding the operations to the list represented as singleton clusters
	#
	# requires that the indicies of `graph` matrix correspond to the indicies of `nodelist`
	# so that distances between nodes can be related to the node objects 
	def build_mst_adjacency_list(graph, nodelist)
		parent = prim(graph) #O(n^2)

		adjacency_list = PriorityQueue.new()
		for i in 1...nodelist.size do #O(n)
			j = parent[i]
			pair = Set[nodelist[i], nodelist[j]]
			adjacency_list.push(pair, graph[i][j]) unless adjacency_list.has_key?(pair) #O(1) (maybe not with this has key check. see build_mst-improvement branch)
		end
		adjacency_list
	end

	def remove_heap_element(heap, obj)
		heap.change_priority(obj, -1)
		heap.delete_min_return_priority
	end

	# this method could be optimized for the case where new_priority is <= old priority
	def replace_heap_element(heap, obj, new_obj, priority, new_priority)
		remove_heap_element(heap, obj)
		heap.push(new_obj, new_priority)
	end

	# naive prims algorithm to find minimum spanning tree
	# O(n^2)
	# @param graph [Array<Array<Float>>]  matrix holding the distance of every pair of nodes in the graph, 
	#                           with a distance of -1 for self-pairings
	# @return [Array<Integer>]  disjoint-set forest where each item id can be used to traverse back up the path travelled 
	def prim(graph)
		n = graph.size
		parent = Array.new(n)
		key = Array.new(n)
		visited = Array.new(n)

		n.times do |i|
			key[i] = Float::MAX
			visited[i] = false
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
		min = Float::MAX
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