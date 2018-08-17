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

# IN ACTIVE DEVELOPMENT, NOT READY FOR USE

needs 'PCR Libs/GradientPcrRepresentation'
needs 'PCR Libs/GradientPcrHelpers'

module GradientPcrBatching
	include GradientPcrHelpers
	include GradientPcrRepresentation

	# Thermocycler information
	CYCLER_COUNT = 4
	ROW_COUNT = 8
	COLUMN_COUNT = 12
	TEMP_RANGE = 10 # degree C allowable gradient temperature range in one thermocycler

	# difference evaluation betweeen 2 extension groups
	# with any less difference than this, extension cluster combination is forced even
	# if it will leave some thermocyclers open
	MANDATORY_EXTENSION_COMBINATION_DIFFERENCE = 3.0 

	# difference evaluation betweeen 2 tanneal groups
	# with any less difference than this, tanneal cluster combination is forced even
	# if it will leave some rows open within a thermocycler
	MANDATORY_TANNEAL_COMBINATION_DIFFERENCE = 0.3 # difference evaluation betweeen 2 groups

	# Given a list of pcr_operations, batches them into CYCLER_COUNT
	# reaction groups, and within each reaction group, batches operations
	# into ROW_COUNT temperature groups, to take advantage of gradient pcr
	#
	# @param [Array<PcrOperation>]  list of operations to find a batching for
	# @param [Array<PcrOperation>]  new (copied) list, where each pcr operation 
	# 						now has batch numbers for both thermocycler and row
	def batch(pcr_operations)

		# Hashmap which will encode groupings of pcr operations
		extension_cluster_to_tanneal_clusters = Hash.new

		extension_clusters = cluster_by_extension_time(pcr_operations) #O(n^2logn)

		return log_clusters(extension_clusters)

		# extension_clusters.each do |extension_cluster|
		# 	tanneal_clusters = cluster_by_annealling_temp(extension_cluster.members)
		# 	extension_cluster_to_tanneal_clusters[extension_cluster] = tanneal_clusters.cluster
		# end

		# result = label_pcr_operations_by_group(pcr_operations, extension_cluster_to_tanneal_clusters)
		# result
	end

	# O(n^2Logn)
	def cluster_by_extension_time(pcr_operations)
		extension_graph = ExtensionClusterGraph.new({
							pcr_operations: 			pcr_operations,
							thermocycler_quantity: 		CYCLER_COUNT,
							thermocycler_rows: 			ROW_COUNT,
							thermocycler_columns: 		COLUMN_COUNT,
							thermocycler_temp_range: 	TEMP_RANGE
						}) #O(n^2)
		extension_graph.checkrep
		while !extension_graph.threshhold_func(MANDATORY_EXTENSION_COMBINATION_DIFFERENCE) #O(n^2logn) for whole loop
			extension_graph.combine_nearest_clusters #O(nlogn)
			extension_graph.checkrep
		end
		extension_graph.cluster_set
	end

	def cluster_by_annealling_temp(pcr_operations)
		tanneal_graph = TannealClusterGraph.new(pcr_operations)
		while !tanneal_graph.threshhold_func(MANDATORY_TANNEAL_COMBINATION_DIFFERENCE)
			tanneal_graph.combine_nearest_clusters
		end
		tanneal_graph.cluster_set
	end

	def label_pcr_operations_by_group(pcr_operations, ext_to_temp_grouping_map)
		batched_pcr_operations = Array.new(pcr_operations.length)
		ext_to_temp_grouping_map.each_with_index do |extension_cluster, tanneal_clusters, ext_idx|
			tanneal_clusters.each_with_index do |tanneal_cluster, temp_idx|
				tanneal_cluster.members.each do |pcr_operation|
					pcr_with_batching_info = pcr_operation.clone()
					pcr_with_batching_info.extension_group = ext_idx
					pcr_with_batching_info.tanneal_group = temp_idx
					batched_pcr_operations << pcr_with_batching_info
				end
			end
		end
		batched_pcr_operations
	end

	def log_clusters(cluster_list)
		result = "#{cluster_list.size} total clusters \n"
		cluster_list.each do |cluster|
			result += "{ " + cluster.to_string + " }\n"
		end
		result
	end
end