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

	# find the range of a pcr operation fields in two clusters
	# if they were to be combined into one
	def combine_range(a,b)
		Max((yeild(a) - yeild(b).abs), (yeild(a) - yeild(b)).abs)
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

	def build_mst_adjacency_list(graph, nodelist)
		parent = prim(graph) #O(n^2)

		adjacency_list = PriorityQueue.new()
		for i in 1...nodelist.size do #O(n)
			j = parent[i]
			pair = [yield(nodelist[i]), yield(nodelist[j])].sort #sorting ensures equality of arrays if same contents
			adjacency_list.push(pair, graph[i][j]) #O(1)
		end
		adjacency_list
	end

	def remove_heap_element(heap, obj)
		heap.change_priority(obj, -1)
		heap.delete_min_return_priority
	end

	# this method could be optimized for the case where new_priority is <= old priority
	def replace_heap_element(heap, obj, new_obj, priority, new_priority)
		heap.remove_heap_element(obj)
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