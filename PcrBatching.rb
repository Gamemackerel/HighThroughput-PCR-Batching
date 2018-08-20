


needs 'PCR Libs/GradientPcrRepresentation'
needs 'PCR Libs/GradientPcrHelpers'


# PcrBatching uses a nearest neighbor chain algorithm to 
# batch pcr operations with their nearest neighbors by extension time until
# CYCLER_COUNT groups remain. 
# These groups are marked by the maximum extension they contain.
# Then, for each thermocycler group, the operations within that group
# are clustered again with a nearest neighbor chain algorithm by temperature, 
# that refuses to make temp groups of larger than COLUMN_COUNT
# until ROW_COUNT temp groups remain.
# each row-grouping is marked with its average temperature.
# range of temp gradient is determined by the highest and lowest
# temp groups. Groups are placed into the rows that have nearest
# temperature to their average, rounding up always.
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
	MANDATORY_EXTENSION_COMBINATION_DIFFERENCE = 30 

	# do not group operations into the same thermocycler if the difference between their extension
	# times is greater than this value
	MAXIMUM_EXTENSION_COMBINATION_DIFFERENCE = 5 * 60

	# difference evaluation betweeen 2 tanneal groups
	# with any less difference than this, tanneal cluster combination is forced even
	# if it will leave some rows open within a thermocycler
	MANDATORY_TANNEAL_COMBINATION_DIFFERENCE = 0.3

	# do not group operations into the same thermocycler row if the difference 
	# between their annealling temperatures is greater than this value
	MAXIMUM_TANNEAL_COMBINATION_DIFFERENCE = 5 * 60


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

		extension_graph = ExtensionClusterGraph.new({
								pcr_operations: 			pcr_operations,
								thermocycler_quantity: 		CYCLER_COUNT,
								thermocycler_rows: 			ROW_COUNT,
								thermocycler_columns: 		COLUMN_COUNT,
								thermocycler_temp_range: 	TEMP_RANGE,
								force_combination_distance: MANDATORY_EXTENSION_COMBINATION_DIFFERENCE,
								prevent_combination_distance: MAXIMUM_EXTENSION_COMBINATION_DIFFERENCE
							}) #O(n^2)

		extension_clusters = extension_graph.perform_clustering #O(n^2logn)
		extension_clusters.each do |extension_cluster|
			tanneal_graph = TannealClusterGraph.new({
									pcr_operations: 			extension_cluster.members,
									thermocycler_rows: 			ROW_COUNT,
									thermocycler_columns: 		COLUMN_COUNT,
									force_combination_distance: MANDATORY_TANNEAL_COMBINATION_DIFFERENCE,
									prevent_combination_distance: MAXIMUM_TANNEAL_COMBINATION_DIFFERENCE
								})

			tanneal_clusters = tanneal_graph.perform_clustering
			extension_cluster_to_tanneal_clusters[extension_cluster] = tanneal_clusters
		end

		extension_cluster_to_tanneal_clusters
	end

	def log_clusters(cluster_list)
		result = "#{cluster_list.size} total clusters \n"
		cluster_list.each do |cluster|
			result += "{ " + cluster.to_string + " }\n"
		end
		result
	end
end