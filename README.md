## High Throughput PCR Batching
Ruby Gem intended for biologists, to make high-throughput PCR easy. Allows optimized batching of as many pcr operations as possible in a set amount of thermocyclers.

Install gem with `gem "pcr_batching", :git => "git://github.com:Gamemackerel/PCR-Batching.git"`

To use: require pcr_batching, then instantiate a PcrBatcher with settings that are appropriate for your lab needs and thermocyclers.

```
    require 'pcr_batching'
    my_batcher = PcrBatcher.new(
                cycler_count: 4,
                row_count: 8,
                column_count: 12,
                temp_range: 17,
                mand_ext_comb_diff: 30.0,
                mand_tanneal_comb_diff: 0.3,
                max_ext_comb_diff: 300.0,
                max_tanneal_comb_diff: 3.0
            )
```

To learn more about these settings, see the [initialize method in PcrBatcher](https://gamemackerel.github.io/PCR-Batching/doc/PcrBatcher.html#initialize-instance_method)

Next, give the PcrBatcher information about each PCR reaction you would like to run using add_pcr_operation. Let's say we need to run 3 reactions at once, two of which have a very similar extension time requirement.

```
    my_batcher.add_pcr_operation(
                extension_time:     59,
                anneal_temp:        69,
                unique_id:          1
            )
    my_batcher.add_pcr_operation(
                extension_time:     60,
                anneal_temp:        74,
                unique_id:          1
            )
    my_batcher.add_pcr_operation(
                extension_time:     600,
                anneal_temp:        70,
                unique_id:          2
            )
```

Finally, we can let the batcher organize these into space efficient groupings. With 4 open thermocyclers and only 3 reactions, this shouldn't be too much trouble. Because the extension times are very similar for pcr operation 1 and 2, we should expect the batcher to place them in the same thermocycler group.

```
    batching_result = my_batcher.batch
```


The result of the batching is a map from extension time group to a list of annealing temperature groups, with the idea that each group of reactions with similar extension time can be run in the same thermocycler, while each subgroup of reactions with similar annealing temperature within that extension time group can be run in the same row of that thermocycler.

If that didn't make sense, let's break down our batching_result to see what this looks like in practice.

```
    key, value = batching_result.first

    key #=> ExtensionCluster<@size: 2>
    key.mean_extension #=> 59.5
    key.members #=> [PcrOperation<@unique_id: 1>, PcrOperation<@unique_id: 2>]
    value #=> [TannealCluster<@size: 1>, TannealCluster<@size: 1>]
    value.first.mean_anneal #=> 69
    value.first.members #=> [PcrOperation<@unique_id: 1>]
```

From the values so far inspected from the hash returned by PcrBatcher, it seems that Pcr 1 and 2 should go in the same thermocycler using an extension time of 59.5 seconds; and within that thermocycler, PcrOperation 1 should be placed on a row close to 69 degrees C.

[View other useful methods in the yard generated documentation](https://gamemackerel.github.io/PCR-Batching/doc).

[Source code available on github](https://github.com/Gamemackerel/PCR-Batching)