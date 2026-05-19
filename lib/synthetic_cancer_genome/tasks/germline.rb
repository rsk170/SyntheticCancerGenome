module SyntheticCancerGenome
  HG19_ORGANISM = Organism.organism_for_build("hg19") || "Hsa/feb2014"
  HG38_ORGANISM = Organism.organism_for_build("hg38") || "Hsa/feb2023"
  GENOMES1000_RELEASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
  GENOMES1000_LOCAL_VCF = File.expand_path("~/.rbbt/share/databases/genomes_1000/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz")
  HG19_TO_HG38_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  HG19_TO_HG38_CHAIN_FILE = File.expand_path("~/.rbbt/share/lift_over/hg19ToHg38.over.chain.gz")
  LIFTOVER_BIN = File.expand_path("~/.rbbt/software/opt/bin/liftOver")

  def self.genomes1000_phase3_vcf_file
    unless File.exist?(GENOMES1000_LOCAL_VCF)
      Open.mkdir File.dirname(GENOMES1000_LOCAL_VCF)
      Open.download GENOMES1000_RELEASE_URL, GENOMES1000_LOCAL_VCF
    end

    GENOMES1000_LOCAL_VCF
  end

  def self.genomes1000_phase3_vcf
    Open.open(genomes1000_phase3_vcf_file)
  end

  def self.hg19_to_hg38_chain_file
    if ! File.exist?(HG19_TO_HG38_CHAIN_FILE) || File.size(HG19_TO_HG38_CHAIN_FILE) < 100_000
      Open.mkdir File.dirname(HG19_TO_HG38_CHAIN_FILE)
      Open.download HG19_TO_HG38_CHAIN_URL, HG19_TO_HG38_CHAIN_FILE
    end

    HG19_TO_HG38_CHAIN_FILE
  end

  input :haploid, :boolean, "Use haploid frequencies (half)", false
  task :genotype_germline_hg19_all_chr => :array do |haploid|
    Open.open_pipe do |sin|
      # The packaged Genomes1000 rsids producer is brittle in this environment.
      # Stream the source VCF directly and extract only the fields needed here.
      vcf_stream = SyntheticCancerGenome.genomes1000_phase3_vcf

      begin
        while line = vcf_stream.gets
          next if line.start_with?("#")

          chr, pos, _id, ref, alt_l, _qual, _filter, info, *_rest = line.chomp.split("\t")
          next if chr.nil? || pos.nil? || ref.nil? || alt_l.nil? || info.nil?

          af_match = info.match(/(?:^|;)EUR_AF=([^;]+)/)
          next if af_match.nil?

          af = af_match[1].to_f
          af = haploid ? af / 2 : af
          next unless rand < af

          pos, alts = Misc.correct_vcf_mutation(pos.to_i, ref, alt_l)
          alt = alts.sample
          next if alt.nil?

          sin.puts [chr, pos, alt] * ":"
        end
      ensure
        vcf_stream.close if vcf_stream.respond_to?(:close) && ! vcf_stream.closed?
      end
    end
  end

  input :chromosome, :string, "Chromosome to choose from", nil
  dep :genotype_germline_hg19_all_chr
  task :genotype_germline_hg19 => :array do |chr|
    chr = chr.sub('chr', '') if chr
    TSV.traverse step(:genotype_germline_hg19_all_chr), :type => :array, :into => :stream do |mutation|
      next if chr && mutation.split(":").first.sub('chr','') != chr
      mutation
    end
  end

  dep :genotype_germline_hg19
  task :genotype_germline_hg38_lf => :array do
    input_count = 0

    TmpFile.with_file do |source_bed|
      TmpFile.with_file do |target_bed|
        TmpFile.with_file do |unmapped_bed|
          Open.write(source_bed) do |bed|
            TSV.traverse step(:genotype_germline_hg19), :type => :array do |position|
              chr, pos, *_rest = position.split(":")
              next if chr.nil? || pos.nil?

              chr = chr.sub(/^chr/, '')
              bed.puts ["chr#{chr}", pos.to_i - 1, pos.to_i, position] * "\t"
              input_count += 1
            end
          end

          chain_file = SyntheticCancerGenome.hg19_to_hg38_chain_file
          CMD.cmd_log("'#{LIFTOVER_BIN}' '#{source_bed}' '#{chain_file}' '#{target_bed}' '#{unmapped_bed}'")
          raise "liftOver produced no mapped positions out of #{input_count} inputs" if input_count > 0 && File.size(target_bed) == 0

          stream = Open.open_pipe do |sin|
            Open.read(target_bed) do |line|
              chr, _start, pos, original = line.chomp.split("\t")
              next if chr.nil? || pos.nil? || original.nil?

              chr = chr.sub(/^chr/, '')
              _old_chr, _old_pos, *rest = original.split(":")
              sin.puts ([chr, pos] + rest) * ":"
            end
          end

          stream
        end
      end
    end
  end

  dep :genotype_germline_hg38_lf
  task :genotype_germline_hg38_lf_chr => :array do 
    chr = recursive_inputs[:chromosome]
    TSV.traverse step(:genotype_germline_hg38_lf), :type => :array, :into => :stream do |line|
      next if ! line.include?(":") || line.split(":").first.include?("_")
      next if chr && line.split(":").first != chr
      line
    end
  end

  dep :genotype_germline_hg38_lf_chr
  dep Sequence, :reference, :positions => :genotype_germline_hg38_lf_chr, :organism => HG38_ORGANISM, :vcf => false, :full_reference_sequence => false
  task :genotype_germline_hg38 => :array do 
    TSV.traverse step(:reference), :into => :stream do |mutation, reference|
      # Make sure we don't take positions that are now non-mutations, as this
      # breaks other tools downstream
      next if mutation.split(":")[2].split(",").include? reference
      mutation
    end
  end

  dep :genotype_germline_hg38, jobname: "father", haploid: true, compute: :produce
  dep :genotype_germline_hg38, jobname: "mother", haploid: true, compute: :produce
  task :diploid_genotype_germline_hg38 => :array do
    Open.open_pipe do |sin|
      TSV.traverse dependencies.first, type: :array, bar: true do |mutation|
        chr, pos, alt = mutation.split(":")
        chr = "chr" + chr unless chr.start_with?("chr")
        chr = "copy-1_" + chr
        sin.puts [chr, pos, alt] * ":"
      end
      TSV.traverse dependencies.last, type: :array, bar: true do |mutation|
        chr, pos, alt = mutation.split(":")
        chr = "chr" + chr unless chr.start_with?("chr")
        chr = "copy-2_" + chr
        sin.puts [chr, pos, alt] * ":"
      end
    end
  end
end
