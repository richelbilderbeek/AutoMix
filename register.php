<?php
include 'db.php';

  // Define post fields into simple variables
  $first_name = $_POST['first_name'];
  $last_name = $_POST['last_name'];
  $email_address = $_POST['email_address'];
  $organisation = $_POST['organisation'];

  /* Let's strip some slashes in case the user entered
any escaped characters. */

  $first_name = stripslashes($first_name);
  $last_name = stripslashes($last_name);
  $email_address = stripslashes($email_address);
  $birth_town = stripslashes($birth_town);
  $username = stripslashes($username);

  /* Do some error checking on the form posted fields */

  if((!$first_name) || (!$last_name) || (!$email_address)){
    $register_ok=0;
  }
  else{	
      /* Everything has passed error checks that we have done.
It's time to create the account! */
      
      // Enter info into the Database.

      $sql = mysql_query("INSERT INTO automix (first_name, last_name, email_address, organisation, signup_date) VALUES('$first_name', '$last_name', '$email_address', '$organisation', now())") or die (mysql_error());

      if(!$sql){
	$register_ok=-2;
	$email_address2="d_hastie@hotmail.com";
	$subject1="AutoMix: Registration failure";   
	$message1="The following registration failed

      first name:    $first_name
      last name:     $last_name
      e-mail:        $email_address
      organisation:  $organisation";

      	mail($email_address2, $subject1, $message1, "From: David Hastie<d_hastie@hotmail.com>\nX-Mailer: PHP/" . phpversion()) or print("Failure E-mail not sent");
          

      } 
      else {
	$register_ok=1;
	$userid = mysql_insert_id();
	// Let's mail the user!
	$subject = "AutoMix";
	$message = "Dear $first_name $last_name,

This is an automated response, please do not reply.

Thank you for downloading the AutoMix software.
	
Your name and email address will be stored in our user database to ensure that we can
send you notification of any significant upgrades to the AutoMix software.

Kind regards
David Hastie";

	
	mail($email_address, $subject, $message, "From: David Hastie<d_hastie@hotmail.com>\nX-Mailer: PHP/" . phpversion()) or print("User notification e-mail not sent. Please e-mail d_hastie@hotmail.com with this message.");

	$email_address2="d_hastie@hotmail.com";
	$subject1="AutoMix: $first_name $lastname has registered to use the software";   
	$message1="New user.";
	
	mail($email_address2, $subject1, $message1, "From: David Hastie<d_hastie@hotmail.com>\nX-Mailer: PHP/" . phpversion()) or print("Author notification e-mail not sent. Please e-mail d_hastie@hotmail.com with this message.");
      }
    }
  
if($register_ok==1){
  include '../php/automix_success.php';
}
else{
  include '../php/automix_failure.php';
}
?>