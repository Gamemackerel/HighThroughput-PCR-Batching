require 'pcr_batching_representation'
require 'pcr_batching_helpers'

# Pcrbatcher uses a nearest neighbor chain algorithm to 
# batch pcr operations with their nearest neighbors by extension time until
# cycler_count groups remain. 
# These groups are marked by the maximum extension they contain.
# Then, for each thermocycler group, the operations within that group
# are clustered again with a nearest neighbor chain algorithm by temperature, 
# that refuses to make temp groups of larger than @column_count
# until @row_count temp groups remain.
# each row-grouping is marked with its average temperature.
# range of temp gradient is determined by the highest and lowest
# temp groups. Groups are placed into the rows that have nearest
# temperature to their average, rounding up always.
#
class PcrBatcher

    # Initialize the important fields needed to batch some pcr operations.
    # @param opts [Hash]  initialization options, will use defaults when option
    #                   is not supplied
    # @option opts [Integer] :cycler_count  number of thermocyclers 
    #                   available for batching
    # @option opts [Integer] :row_count  number of rows per thermocycler
    # @option opts [Integer] :column_count  number of columns per thermocycler
    # @option opts [Float] :temp_range  degree C allowable gradient temperature
    #                   range in one thermocycler
    # @option opts [Float] :mand_ext_comb_diff  the mandatory extension
    #                   combination difference at which pcr operations will be
    #                   placed into the same thermocycler and use the same
    #                   extension time, even if it will leave
    #                   some thermocyclers open
    # @option opts [Float] :max_ext_comb_diff  the extension difference at
    #                   which pcr operations will be disallowed from being 
    #                   placed into the same thermocycler and using the same 
    #                   extension time, even if there will be
    #                   more pcr batches than thermocyclers
    # @option opts [Float] :mand_tanneal_comb_diff  the mandatory tanneal
    #                   combination difference at which pcr operations will be
    #                   placed into the same thermocycler row, even that means
    #                   some rows of the thermocylcer will be empty
    # @option opts [Float] :max_tanneal_comb_diff  the tanneal difference at
    #                   which pcr operations will be disallowed from being
    #                   grouped into the same thermocycler row
    def initialize opts = {}
        opts = defaults.merge opts
        @cycler_count           = opts[:cycler_count]
        @row_count              = opts[:row_count]
        @column_count           = opts[:column_count]
        @temp_range             = opts[:temp_range]
        @mand_ext_comb_diff     = opts[:mand_ext_comb_diff]
        @mand_tanneal_comb_diff = opts[:mand_tanneal_comb_diff]
        @max_ext_comb_diff      = opts[:max_ext_comb_diff]
        @max_tanneal_comb_diff  = opts[:max_tanneal_comb_diff]
        @pcr_operations         = []
    end

    # Batching settings that work well for the UW BIOFAB PCR workflow
    # which uses 4 X thermocyclers
    def defaults
        {
            cycler_count: 4,
            row_count: 8,
            column_count: 12,
            temp_range: 17,
            mand_ext_comb_diff: 30.0,
            mand_tanneal_comb_diff: 0.3,
            max_ext_comb_diff: 300.0,
            max_tanneal_comb_diff: 3.0,
        }
    end

    # Make a brand new pcr operation representation and add it to the list
    # of pcr operations to be batched.
    #
    # @param opts [Hash]  pcr operation definition options
    # @option opts [Float] :extension_time  extension time for this
    #               pcr operation
    # @option opts [Float] :anneal_temp  the annealing temperature for
    #               this pcr operation
    # @option opts [Integer] :extension_group  a group id for this operation,
    #               shared with other pcr operations who could be run together
    #               in a reaction with the same extension time
    # @option opts [Integer] :tanneal_group  a group id for this operation,
    #               shared with other pcr operations who could be run together
    #                in a reaction with the same annealing temperature
    # @option opts [Integer] :unique_id  identifier to track a pcr operation
    #                through its being batched
    def add_pcr_operation opts = {}
        @pcr_operations << PcrOperation.new(
            extension_time:  opts[:extension_time],
            anneal_temp:     opts[:anneal_temp],
            unique_id:       opts[:unique_id],
        )
    end

    # Batches pcr_operations into cycler_count
    # reaction groups, and within each reaction group, batches operations
    # into @row_count temperature groups.
    # 
    # @param use_checkrep [Boolean]  whether or not to check representation invariants
    #               during clustering; very slow and only necessary while testing
    # @return [Hash<ExtensionCluster, Set<TannealCluster>>]  a mapping from 
    #               a group of pcr operations with similar extension time to the set of  
    #               sub-groups of that group which have similar anneal temperature. 
    def batch(use_checkrep = false)
        extension_cluster_to_tanneal_clusters = Hash.new

        extension_graph = ExtensionClusterGraph.new(
                    pcr_operations:               @pcr_operations,
                    thermocycler_quantity:        @cycler_count,
                    thermocycler_rows:            @row_count,
                    thermocycler_columns:         @column_count,
                    thermocycler_temp_range:      @temp_range,
                    force_combination_distance:   @mand_ext_comb_diff,
                    prevent_combination_distance: @max_ext_comb_diff
                ) # O(n^2)

        extension_clusters = extension_graph.perform_clustering use_checkrep # O(n^2)
        extension_clusters.each do |extension_cluster| # m clusters =>  O(m) 
            tanneal_graph = TannealClusterGraph.new(  # q pcr operations per cluster => O(q^2) 
                        pcr_operations:               extension_cluster.members,
                        thermocycler_rows:            @row_count,
                        thermocycler_columns:         @column_count,
                        force_combination_distance:   @mand_tanneal_comb_diff,
                        prevent_combination_distance: @max_tanneal_comb_diff
                    )

            tanneal_clusters = tanneal_graph.perform_clustering use_checkrep #O(q^2)
            extension_cluster_to_tanneal_clusters[extension_cluster] = tanneal_clusters
        end # m * q = n => O(2 * n^2) 

        extension_cluster_to_tanneal_clusters
    end # O(4 * n^2)
end