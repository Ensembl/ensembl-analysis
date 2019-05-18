# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use warnings ;
use strict;

package Bio::EnsEMBL::Analysis::Tools::IMGT::SeqIO::imgt_embl;

use base qw(Bio::SeqIO::embl);

sub _initialize {
  my($self,@args) = @_;
  
  $self->SUPER::_initialize(@args);

  $self->sequence_factory(new Bio::Seq::SeqFactory
                          (-verbose => $self->verbose(),
                           -type => 'Bio::Seq::RichSeqIMGT'));

}



sub next_seq {
  my ($self,@args) = @_;
  my ($pseq,$c,$line,$name,$dataclass,$desc,$acc,$seqc,$mol,$div,
      $date, $comment, @date_arr);
  
  my ($annotation, %params, @features) =
      new Bio::Annotation::Collection;
  
  $line = $self->_readline;
  # This needs to be before the first eof() test
  
  if( !defined $line ) {
    return; # no throws - end of file
  }
  
  if( $line =~ /^\s+$/ ) {
    while( defined ($line = $self->_readline) ) {
      $line =~/^\S/ && last;
    }
    # return without error if the whole next sequence was just a single
    # blank line and then eof
    return unless $line;
  }
  
  # no ID as 1st non-blank line, need short circuit and exit routine
  $self->throw("EMBL stream with no ID. Not embl in my book")
      unless $line =~ /^ID\s+\S+/;
  
  # At this point we are sure that $line contains an ID header line
  my $alphabet;
  # IMGT header : ID   A00673 IMGT/LIGM annotation : keyword level; DNA; SYN; 45 BP.
  if ( $line =~ /^ID\s+(\S+)\s+\S+\s+\S+\s+:\s+([^;]+);\s+([^;]+);\s+([^;]+);\s+\d+\s+BP\.$/) {   
    ($name, $dataclass, $mol, $div) = ($1, $2, $3, $4);

    if (defined $mol ) {
      if ($mol =~ /DNA/) {
        $alphabet='dna';
      }
      elsif ($mol =~ /RNA/) {
        $alphabet='rna';
      }
      elsif ($mol =~ /AA/) {
        $alphabet='protein';
      }
    }
  }
  
  unless( defined $name && length($name) ) {
    $name = "unknown_id";
  }
  
  # $self->warn("not parsing upper annotation in EMBL file yet!");
  my $buffer = $line;
  local $_;
  BEFORE_FEATURE_TABLE: until( !defined $buffer ) {
    $_ = $buffer;
    # Exit at start of Feature table
    if( /^(F[HT]|SQ)/ ) {
      $self->_pushback($_) if( $1 eq 'SQ' );
      last;
    }
    # Description line(s)
    if (/^DE\s+(\S.*\S)/) {
      $desc .= $desc ? " $1" : $1;
    }
    
    #accession number
    if( /^AC\s+(.*)?/ ) {
      my @accs = split(/[; ]+/, $1); # allow space in addition
      $params{'-accession_number'} = shift @accs
          unless defined $params{'-accession_number'};
      push @{$params{'-secondary_accessions'}}, @accs;
    }
    
    #version number
    if( /^SV\s+\S+\.(\d+);?/ ) {
      my $sv = $1;
      #$sv =~ s/\;//;
      $params{'-seq_version'} = $sv;
      $params{'-version'} = $sv;
    }
    
    #date (NOTE: takes last date line)
    if( /^DT\s+(.+)$/ ) {
      my $line = $1;
      my ($date, $version) = split(' ', $line, 2);
      $date =~ tr/,//d; # remove comma if new version
      if ($version =~ /\(Rel\. (\d+), Created\)/xms ) {
        my $release = Bio::Annotation::SimpleValue->
            new(
                -tagname    => 'creation_release',
                -value      => $1
                );
        $annotation->add_Annotation($release);
      } elsif ($version =~ /\(Rel\. (\d+), Last updated, Version (\d+)\)/xms ) {
        my $release = Bio::Annotation::SimpleValue->
            new(
                -tagname    => 'update_release',
                -value      => $1
                );
        $annotation->add_Annotation($release);
        
        my $update = Bio::Annotation::SimpleValue->
            new(
                -tagname    => 'update_version',
                -value      => $2
                );
        $annotation->add_Annotation($update);
      }
      push @{$params{'-dates'}}, $date;
    }
    
    #keywords
    if( /^KW   (.*)\S*$/ ) {
      my @kw = split(/\s*\;\s*/,$1);
      push @{$params{'-keywords'}}, @kw;
    }
    
    # Organism name and phylogenetic information
    elsif (/^O[SC]/) {
      # pass the accession number so we can give an informative throw message if necessary
      my $species = $self->_read_EMBL_Species(\$buffer, $params{'-accession_number'});
      $params{'-species'}= $species;
    }
    
    # References
    elsif (/^R/) {
      my @refs = $self->_read_EMBL_References(\$buffer);
      foreach my $ref ( @refs ) {
        $annotation->add_Annotation('reference',$ref);
      }
    }
    
    # DB Xrefs
    elsif (/^DR/) {
      my @links = $self->_read_EMBL_DBLink(\$buffer);
      foreach my $dblink ( @links ) {
        $annotation->add_Annotation('dblink',$dblink);
      }
    }
    
    # Comments
    #elsif (/^CC\s+(.*)/) {
    #  $comment .= $1;
    #  $comment .= " ";
    #  while (defined ($_ = $self->_readline) ) {
    #    if (/^CC\s+(.*)/) {
    #      $comment .= $1;
    #      $comment .= " ";
    #    }
    #    else {
    #      last;
    #    }
    #  }
    #  my $commobj = Bio::Annotation::Comment->new();
    #  $commobj->text($comment);
    #  $annotation->add_Annotation('comment',$commobj);
    #  $comment = "";
    #}

    # Get next line.
    $buffer = $self->_readline;
  }
  
  while( defined ($_ = $self->_readline) ) {
    /^FT\s{3}\w/ && last;
    /^SQ / && last;
    /^CO / && last;
  }
  $buffer = $_;
  
  if (defined($buffer) && $buffer =~ /^FT\s+\S+/) {
    until( !defined ($buffer) ) {
      
      my $ftunit = $self->_read_FTHelper_EMBL(\$buffer);
      next if not defined $ftunit;

      # process ftunit
      my $feat =
          $ftunit->_generic_seqfeature($self->location_factory(), $name);
      
      # add taxon_id from source if available
      if($params{'-species'} && ($feat->primary_tag eq 'source')
         && $feat->has_tag('db_xref')
         && (! $params{'-species'}->ncbi_taxid())) {
        foreach my $tagval ($feat->get_tag_values('db_xref')) {
          if(index($tagval,"taxon:") == 0) {
            $params{'-species'}->ncbi_taxid(substr($tagval,6));
            last;
          }
        }
      }
      
      # add feature to list of features
      push(@features, $feat);
      
      if( $buffer !~ /^FT/ ) {
        last;
      }
    }
  }
  # skip comments
  while( defined ($buffer) && $buffer =~ /^XX/ ) {
    $buffer = $self->_readline();
  }
  
  if( $buffer !~ /^SQ/  ) {
    while( defined ($_ = $self->_readline) ) {
      /^SQ/ && last;
    }
  }
  $seqc = "";
  while( defined ($_ = $self->_readline) ) {
    m{^//} && last;
    $_ = uc($_);
    s/[^A-Za-z]//g;
    $seqc .= $_;
  }
  my $seq = $self->sequence_factory->create
      (-verbose => $self->verbose(),
       -division => $div,
       -data_class => $dataclass,
       -seq => $seqc,
       -desc => $desc,
       -display_id => $name,
       -primary_id => $name,
       -annotation => $annotation,
       -molecule => $mol,
       -alphabet => $alphabet,
       -features => \@features,
       %params);
  return $seq;
}



=head2 _read_EMBL_Species

 Title   : _read_EMBL_Species
 Usage   :
 Function: Reads the EMBL Organism species and classification
           lines.
 Example :
 Returns : A Bio::Species object
 Args    : a reference to the current line buffer, accession number
 Notes   : Over-ridden from Bio/SeqIO/embl.pm to cope with "eccentric"
           OS and OC lines in IMGT flat file
=cut

sub _read_EMBL_Species {
  my( $self, $buffer, $acc ) = @_;
  my $org;
  
  $_ = $$buffer;
  my( $sub_species, $species, $genus, $common, $sci_name, $class_lines );
  while (defined( $_ ||= $self->_readline )) {
    if (/^OS\s+(.+)/) {
      $sci_name .= ($sci_name) ? ' '.$1 : $1;
    }
    elsif (s/^OC\s+(.+)$//) {
      $class_lines .= $1;
    }
    elsif (/^OG\s+(.*)/) {
      $org = $1;
    }
    else {
      last;
    }
    
    $_ = undef; # Empty $_ to trigger read of next line
  }

  $$buffer = $_;
  
  $sci_name || return;

  # IMGT sometimes has "unclassified" for both species and taxon
  # which confuses Bio::Taxon. Replace it with unidentified in species
  $sci_name =~ s/unclassified/unidentified/i;

  # sometimes things have common name in brackets, like
  # Schizosaccharomyces pombe (fission yeast), so get rid of the common
  # name bit. Probably dangerous if real scientific species name ends in
  # bracketed bit.
  if ($sci_name =~ /^(.+)\s+\((.+)\)$/) {
    $sci_name = $1;
    $common  = $2;
  }

  if ($sci_name =~ /^(\S+)\s+(\S+)$/) {
    ($genus, $species) = ($1, $2);
  } elsif ($sci_name =~ /(\S+)\s+(\S+)\s+(\S+)/) {
    ($genus, $species, $sub_species) = ($1, $2, $3);
    $sci_name = $genus . " " . $species;
  }

  # Convert data in classification lines into classification array.
  # only split on ';' or '.' so that classification that is 2 or more words
  # will still get matched, use map() to remove trailing/leading/intervening
  # spaces
  my @class = map { s/^\s+//; s/\s+$//; s/\s{2,}/ /g; $_; } split /[;\.]+/, $class_lines;
  @class = grep { /\S/ } @class;
  
  # Bio::Species array needs array in Species -> Kingdom direction
  unless ($class[-1] eq $species) {
    push(@class, $sci_name);
  }
  @class = reverse @class;
  
  my %names;
  foreach my $i (0..$#class) {
    my $name = $class[$i];
    $names{$name}++;
    if ($names{$name} > 1 && $name ne $class[$i - 1]) {
      $self->throw("$acc seems to have an invalid species classification.");
    }
  }
  
  my $make = Bio::Species->new();
  $make->scientific_name($sci_name);
  $make->classification(@class);
  $make->genus($genus) if $genus;
  $make->species($species) if $species;
  $make->sub_species($sub_species) if $sub_species;
  $make->common_name($common) if $common;
  $make->organelle($org) if $org;
  return $make;
}

=head2 _read_FTHelper_EMBL

 Title   : _read_FTHelper_EMBL
 Usage   : _read_FTHelper_EMBL($buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::SeqIO::FTHelper object
 Args    : filehandle and reference to a scalar
 Notes   : overridden from Bio::SeqIO::embl
           to cope with blank FT lines in IMGT flatfile
=cut

sub _read_FTHelper_EMBL {
  my ($self,$buffer) = @_;
  
  my ($key,   # The key of the feature
      $loc,   # The location line from the feature
      @qual,  # An arrray of lines making up the qualifiers
      );
  
  if ($$buffer =~ /^FT\s{3}(\S+)\s*(\S*)$/ ) {
    $key = $1; 
    $loc = $2;
    
    if (not $loc) {
      # in some cases, there is no whitespace between
      # the key and the location. Try to separate them.
      if ($key =~ /^(\S+?)([\>\<\.\d]+)/) {
        $key = $1;
        $loc = $2;
      }
    }

    # Read all the lines up to the next feature
    while ( defined($_ = $self->_readline) ) {
      if (/^FT(\s+)(.+?)\s*$/) {
        # Lines inside features are preceeded by 19 spaces
        # A new feature is preceeded by 3 spaces
        my ($spacer, $rest) = ($1, $2);
        next if $rest !~ /\S/;

        if (length($spacer) > 4) {
          # Add to qualifiers if we're in the qualifiers
          if (@qual) {
            push(@qual, $rest);
          }
          # Start the qualifier list if it's the first qualifier
          elsif (substr($rest, 0, 1) eq '/') {
            @qual = ($rest);
          }
          # We're still in the location line, so append to location
          else {
            $loc .= $rest;
          }
        } else {
          # We've reached the start of the next feature
          last;
        }
      } else {
        # We're at the end of the feature table
        last;
      }
    }    
  } else {
    # No feature key
    return;
  }
  
  # Put the first line of the next feature into the buffer
  $$buffer = $_;

  return if not $key or not $loc;
  
  # Make the new FTHelper object
  my $out = new Bio::SeqIO::FTHelper();
  $out->verbose($self->verbose());
  $out->key($key);
  $out->loc($loc);
  
  # Now parse and add any qualifiers.  (@qual is kept
    # intact to provide informative error messages.)
  QUAL: for (my $i = 0; $i < @qual; $i++) {
    $_ = $qual[$i];
    my( $qualifier, $value ) = m{^/([^=]+)(?:=(.+))?}
    or $self->throw("Can't see new qualifier in: $_\nfrom:\n"
                    . join('', map "$_\n", @qual));
    if (defined $value) {
      # Do we have a quoted value?
      if (substr($value, 0, 1) eq '"') {
        # Keep adding to value until we find the trailing quote
        # and the quotes are balanced
      QUOTES: while ($value !~ /\"$/ or $value =~ tr/"/"/ % 2) { #"
        $i++;
        my $next = $qual[$i];
        if (!defined($next)) {
          $self->warn("Unbalanced quote in:\n".join("\n", @qual).
                      "\nAdding quote to close...".
                      "Check sequence quality!");
          $value .= '"';
          last QUOTES;
        }
        
        # Join to value with space if value or next line contains a space
        $value .= (grep /\s/, ($value, $next)) ? " $next" : $next;
      }
        # Trim leading and trailing quotes
        $value =~ s/^"|"$//g;
        # Undouble internal quotes
        $value =~ s/""/"/g; #"
      }
    } else {
      $value = '_no_value';
    }
    
    # Store the qualifier
    $out->field->{$qualifier} ||= [];
    push(@{$out->field->{$qualifier}},$value);
  }
  
  return $out;
}


1;

