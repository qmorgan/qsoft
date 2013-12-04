  #!/bin/sh

  #Backup script for nimrodel Macbook Pro to MacBackup (Cavalry 1 TB harddrive)  

  # To use rsyncx:
  # RSYNC=/usr/local/bin/rsync --eahfs --showtogo 
  # To use built-in rsync (OS X 10.4 and later):
  RSYNC=/usr/bin/rsync 
  
  # sudo runs the backup as root
  # --eahfs enables HFS+ mode
  # -a turns on archive mode (recursive copy + retain attributes)
  # -x don't cross device boundaries (ignore mounted volumes)
  # -S handle sparse files efficiently
  # --showtogo shows the number of files left to process
  # --delete deletes any files that have been deleted locally
  # $* expands to any extra command line options you may give

  sudo $RSYNC -a -x -S -E --delete --progress \
    --exclude-from nimrodel_backup_excludes.txt $* / /Volumes/MacBackup/

  # make the backup bootable - comment this out if needed

  sudo bless -folder /Volumes/MacBackup/System/Library/CoreServices
