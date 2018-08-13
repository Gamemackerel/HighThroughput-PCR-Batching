module GradientPcrHelpers
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
		Max(abs(yeild(a) - yeild(b)), abs(yeild(a) - yeild(b))
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

	def change_priority(heap, obj, new_priority, old_priority)
		if new_priority < old_priority
			heap.decrease_priority(obj, new_priority) # O(1) with good priority list, or O(Logn) with naive one. Fastcontainers is naive, but it is also 5x as fast as alternatives for push and pop
		elsif new_priority > old_priority
			heap.increase_priority(obj, new_priority) # O(Logn)
		end
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