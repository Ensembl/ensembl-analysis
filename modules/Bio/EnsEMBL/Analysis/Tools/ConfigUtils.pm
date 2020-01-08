# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::ConfigUtils;

use warnings ;
use strict;
use FindBin ;
use lib "$FindBin::Bin" ;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
#use Module::Load;
use Data::Dumper ;
use Bio::EnsEMBL::Analysis::Tools::Stashes qw( package_stash ) ;
use PPI;

use vars qw(@ISA) ;
@ISA = qw();

sub new {
  my ($class,@args) = @_;
  my $self = bless {},$class;

  my ($config_module, $is_example) = rearrange (['CONFIG', 'IS_EXAMPLE'], @args);

  $self->config_module($config_module, $is_example) if $config_module;
  $self->{_doc}= PPI::Document->new($self->file($is_example));
  return $self;
}

=head2 write_config

  Example    : $self->write_config()
  Description: writes new config file to disk
               package statement, pod and import-method are read out of old config file
               and transferred to new config file the comments in the Config hash section
               will be lost

=cut


sub write_config {
   my ($self) = @_ ;

   print "backing up old config file if it exists...\n";

   $self->backup_existing_config() if ( -e $self->file ) ;

   my @content;
   push @content, "\n\n", @{$self->get_Pod},"\n\n", @{$self->get_package_statement};
   push @content, "\n\n" , $self->print;
   push @content, "\n\n",$self->get_subroutine("import"),"\n\n1;";

   my $f = $self->file ;
   open (I,">$f") || die "Can't access file ".$self->file()."\n";
     print I join("\n" ,@content) ;
   close(I);
}

=head2 backup_existing_config

  Example    : $self->backup_existing_config()
  Description: tests if the original config file exists, if the file exists it
               copies it to a new file and extends the filename by ".bak.X"

  Returntype : true/false

=cut

sub backup_existing_config {
   my ($self) = @_ ;
   my $backuped_file = $self->file;
   my $count = 0 ;
   while ( -e $backuped_file ) {
     #print "file exists : $backuped_file \n" ;
     $backuped_file = $self->file.".bak.$count" ;
     $count++;
   }
   my $old_file = $self->file;
   #my $cmd = "cp " . $self->file ." $backuped_file";
   #`$cmd`;
   open(OLD,"$old_file") || throw(" could not read old config $old_file\n");
     my @conf= <OLD> ;
   close(OLD) ;

   my $backup_time = `date +%a","%d" "%B" "%T`;

   unshift (@conf, "\#\n\# this file was backuped on $backup_time\#\n\n");

   open(NEW,">$backuped_file") || throw "Could not write to file $backuped_file\n" ;
     print NEW join("", @conf ) ;
   close(NEW) ;

   print "file $old_file backuped + modified to :\n$backuped_file\n" ;
}


sub get_config_hash {
   my ($self) = @_ ;

   my $rs ;
   my $v = $self->varname ;
   for my $s ( @{$self->get_all_statements} ) {
     if ($s->first_token=~m/$v/ && $s->first_token->isa("PPI::Token::Symbol")) {
        # config_hash found : %Config %Databases ...
        throw("statement appears twice in config - error!\n") if $rs ;
        $rs=$s;
     }
   }
  return $rs ;
}


sub get_all_statements{
   my ($self) = @_ ;
   my @rs;
   for my $s (@{$self->{_doc}->find("PPI::Statement")}) {
     push @rs,$s unless ref($s)=~m/PPI::Statement::/;
   }
   return \@rs ;
}

sub get_X{
   my ($self) = @_ ;
   #return $self->{_doc}->find("PPI::Structure::Constructor");
   return $self->{_doc}->find("PPI::Node");
   #return $self->{_doc}->find("PPI::Element");
   #return $self->{_doc}->find("PPI::Structure::List");
}

sub get_Block{
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Structure::List");
}


sub get_all_symbols {
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Token::Symbol");
}

sub get_Comments{
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Token::Comment");
}

sub get_Variables {
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Statement::Variable");
}

sub get_Pod {
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Token::Pod");
}

sub get_package_statement {
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Statement::Package");
}

sub get_include_statements {
   my ($self) = @_ ;
   return $self->{_doc}->find("PPI::Statement::Include");
}


# hands back subroutine with $name out of config file

sub get_subroutine {
   my ($self,$name) = @_ ;

   my @subs = @{$self->{_doc}->find("PPI::Statement::Sub")};
   for my $s ( @subs ) {
     return $s if ($s->name eq $name);
   }
   return undef;
}



=head2 append_href

  Arg [1]    : Hashref
  Arg [2]    : opt. String
  Example    :
  Description: appends a new config value (href, aref) to the configuration hash
               if Arg[2] is given the new key will be placed in this level / section
  Returntype :
  Exceptions : none

=cut

sub append_href{
  my ($self, $href_to_append, $section, $logic_name) = @_ ;

  warning("$logic_name already exists in config file") if ($self->exists_in_config($logic_name));

  # now walk trough hash and replace at the !!correct!! level
  # but how to we get the correct level ? ? ?

  my %new_config = %{ $self->configuration }  ;

  $new_config{$section}{$logic_name} = $href_to_append ;
  $self->configuration(\%new_config) ;
}


sub get_config_by_name {
   my  ($self, $name, $config_href)  = @_ ;

   # this method can either be called on an object itself
   #  $self->get_config_by_name as well as
   #  in a non oo style for recuriveness
   #  get_config_by_name( $obj, $name )

   my %tmp ;
   unless ( $config_href ) {
     %tmp = %{$self->configuration} ;
   } else {
     %tmp = %$config_href ;
   }

   my $result ;
   KEY: foreach my $key ( keys %tmp ) {
      next KEY unless $tmp{$key};

      if ( $key eq $name ) {
         $self->result($tmp{$key}) ;
         $result = $tmp{$key} ;
         return $result;
      }
      if (ref($tmp{$key}) =~m/HASH/  ){
         $self->get_config_by_name($name, $tmp{$key});
      }
   }
  return $result ;
}


sub result {
  my ($self, $r) = @_ ;
  $self->{_result} = $r if $r ;
  return $self->{_result};
}


sub get_all_keys {
  my ($self) = @_ ;
  my @keys ;
  my %tmp = %{$self->configuration};

  my @ak = @{print_keys (\%tmp)} ;
  my %hc ;
  @hc{@ak} =1 ;
  return [keys %hc] ;
}

# returns list of all keys which hold an anonymous hash in config


sub print_keys {
   my ($t) = @_ ;
   my %tmp= %{$t};
   my @ak ;
   KEY: foreach my $key ( %tmp ) {
      next KEY unless $key ;
      if (ref($tmp{$key}) =~m/HASH/  ){
         push @ak, @{print_keys($tmp{$key})};
      }
     push @ak , $key   if (ref($tmp{$key}) =~m/HASH/) ;
  }
  return \@ak;
}


sub get_config_section {
  my ($self,$section) = @_ ;

  my %tmp = %{$self->configuration} ;
  foreach my $k ( keys %tmp ) {
    print " get_config_section : $k\tval $tmp{$k}\n" ;
    if ($k eq $section){
       print "found $k\n" ;
       return $tmp{$k} ;
    }
  }
}


sub exists_in_config{
   my ($self, $name) = @_;
   $self->check();

   my @to_check = @{$self->get_all_keys() } ;
   my %hc ;
   @hc{@to_check} = 1 ;
   return 1 if ( exists $hc{$name} ) ;
   return 0 ;
}


sub check {
  my ($self)  = @_ ;
  throw ( " you need to hand over a Bio::EnsEMBL::Analysis::Tools::ConfigUtils object")
  unless ( $self->isa("Bio::EnsEMBL::Analysis::Tools::ConfigUtils")) ;

}

=head2 config_module

  Arg [1]    : String describing Class name of config
  Example    : $config->config_module("Bio::EnsEMBL::Analysis::Config::Exonerate2Genes");
  Description: getter/setter for config name
  Returntype : string
  Exceptions : none

=cut

sub config_module{
  my ($self, $c, $is_example) =@_ ;
  if ( $c ) {
    $self->{_config} = $c;
    # setenv PERL5LIB /nfs/acari/jhv/lib:${PERL5LIB}
    # load $c;
    #
    # we could use self->file here
    my $d = $c;
    $c .= '.pm';
    $c .= '.example' if ($is_example);
    $c=~ s{::}{/}g;
    require $c;

    my ($config_href, $varname ) = @{package_stash("$d")};
    $self->configuration($config_href) ;
    $self->varname($varname) ;
  }
  return $self->{_config} ;
}




sub print {
  my ($self) = @_ ;

  my $var_name = $self->varname;

  my $d = Data::Dumper->new([$self->configuration],[( "*var_name" )]);
  # $d->Purity(1);
  # $d->Terse(1);
  $d->Indent(1);
  $d->Deepcopy(1);
  my $rep = $d->Dump() ;
  $rep=~s/\%var_name/\%$var_name/;
  return $rep ;
}




=head2 varname

  Arg [1]    : String
  Example    : $config->varname($string)
  Description: getter/setter for hash-name of the variable used in the config file
  Returntype : String
  Exceptions : none

=cut

sub varname{
  my ($self, $c) =@_ ;
  $self->{_varname} = $c if $c ;
  return $self->{_varname} ;
}




=head2 configuration

  Arg [1]    : Hashref
  Example    : $config->configuration(\%config)
  Description: getter/setter for the configuration hash exported by config module
  Returntype : Hashref
  Exceptions : none

=cut

sub configuration {
  my ($self, $c) =@_ ;
  $self->{_configuration} = $c if $c ;
  return $self->{_configuration} ;
}



=head2 file

  Arg [1]    : none
  Example    : $config->file; print $rf->name();
  Description: Returns absolute path to config file
  Returntype : string
  Exceptions : none

=cut


sub file {
  my ( $self, $is_example ) = @_ ;
  (my $path = $self->config_module) =~s/\:\:/\//g;
  $path .= '.pm';
  $path .= '.example' if ($is_example);
  return $INC{$path} ;
}




1;
