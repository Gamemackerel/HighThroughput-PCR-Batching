Gem::Specification.new do |s|
  s.name        = 'pcr_batching'
  s.version     = '0.0.0'
  s.date        = '2018-12-22'
  s.summary     = "Smart Gradient Polymerase Chain Reaction batching"
  s.description = 'A clustering algorithm useful to laboratory biologists for increasing the efficiency of a PCR workflow. Using this batching-aid, many seperate PCR reactions can be performed at once in one or more "Gradient PCR" enabled thermocyclers, by grouping the reactions into coherent batches that use a shared extension time and anneal temperature. Each row of each thermocycler can accomodate a single reaction group. For the source code, refer to `https://github.com/Gamemackerel/PCR-Batching`.'
  s.authors     = ["Abraham Miller"]
  s.email       = 'abemill@uw.eu'
  require 'rake'
  s.files = FileList[ 'lib/*',
                      'doc/*',
                      '*',
                      'test/*'].to_a
  s.homepage    =
    'http://rubygems.org/gems/pcr-batching'
  s.license       = 'MIT'
end